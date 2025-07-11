import json
import numpy as np
import argparse
from string import ascii_uppercase

class AlphaIndex:
    """ ascii_uppercase can run out of letters for large protein complexes, such as the proteasome.
        Then, two-letter combinations are generated.
    """
    def __init__(self):
        self.letters = ascii_uppercase
        self.base = len(self.letters)

    def __getitem__(self, index):
        if index < 0:
            raise IndexError("Negative indexing not supported.")
        result = ""
        while True:
            index, rem = divmod(index, self.base)
            result = self.letters[rem] + result
            if index == 0:
                break
            index -= 1  # Adjust for zero-based indexing
        return result

alpha = AlphaIndex()

def sequence_template(seq_id, seq, use_template=True, use_msa=True):
    """ Print a single seuqence with a given id.
    """

    seq_rec = {
      "protein": {
        "id": seq_id,
        "sequence": seq,
      }
    }

    # Set explicitly to empty string if not used
    if not use_template:
        seq_rec["protein"]["templates"] = ""

    if not use_msa:
        seq_rec["protein"]["unpairedMsa"] = ""
        seq_rec["protein"]["pairedMsa"] = ""

    return seq_rec


def linker_template(seq_id, linker_name="TME"):
    """ Create a linker ligand, defaulting to a propane (TME)
    """
    linker = {
      "ligand": {
        "id": seq_id,
        "ccdCodes": [linker_name]
      }
    }

    return linker


def bondedAtomPairs(cter_aa, mid_aa, link_cter, link_mid, linker_A="C1", linker_B="C3"):
    """
    """
    return [
      [[mid_aa ,  link_mid,  "CB"], ["L"+cter_aa, 1, linker_A]],
      [[cter_aa, link_cter, "OXT"], ["L"+cter_aa, 1, linker_B]]
    ]


class UbiquitinChain(object):

    # WT ubiquitin sequence
    ub_seq = "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"

    def __init__(self, fname, name="polyub", extra_seq=None, use_template=True, use_msa=True):

        # The json dictionary to populate.
        self.out_json = {
            "name": name,
            "modelSeeds": [1,],
            "dialect": "alphafold3",
            "version": 2,
            "sequences": [],
        }
        self.use_template = use_template
        self.use_msa = use_msa

        # Load and store the connectivity
        con = np.loadtxt(fname, dtype=str, ndmin=2)
        self.connections = con
        self.extra_seq = extra_seq

        # Get the unique Ub chains
        ub_ids = np.unique(self.connections[:,:2].flatten())
        # The second, '-' character means no linkage!
        mask = ub_ids == '-'
        self.ub_ids = np.sort(ub_ids[~mask])
        self.n_links = len(self.connections) - mask.sum()

        # Separate the two kinds of linkages
        mask = con[:,2] == '1'
        self.m1_links = con[mask]
        proper_links = con[~mask]
        self.proper_links = proper_links[~(proper_links[:,1] == '-')]

        # Treat the m1 links
        self.connect_m1_chains()

        # Mutate the anchors for the linkages
        self.mutate_anchor()

        # Set up the sequence list
        self.create_sequence_list()

        # Create the bonds
        self.create_bonds()

        with open(f'{name}.json', 'w', encoding='utf-8') as f:
            json.dump(self.out_json, f, indent=4, ensure_ascii=False)


    def connect_m1_chains(self):
        # Create list with single chain IDs, to be combined
        chain_order = self.ub_ids.astype(f"<U{len(self.ub_ids)}")

        # Introduce all the M1 linkages
        for chB, chA, _ in self.m1_links:
            # Look up the current entry with chains A and B
            # NOTE: it's an error if there's no such element!
            for i,elem in enumerate(chain_order):
                if chA in elem:
                    idA  = i
                    seqA = elem
                if chB in elem:
                    idB  = i
                    seqB = elem

            assert seqA != seqB, "Error! Trying to link an already linked chain!"

            chain_order[idA] = seqB + seqA
            chain_order = np.delete (chain_order, idB)

        # Store indices in a pointer and create the sequences
        seq_ptr = {}
        seq_tab = []
        for i, line in enumerate(chain_order):
            for j, elem in enumerate(line):
                seq_ptr[elem] = [i,j*len(self.ub_seq)]
            seq_tab.append(self.ub_seq * len(line))

        self.seq_tab = seq_tab
        self.seq_ptr = seq_ptr


    def mutate_anchor(self):
        """ Modify the Ub sequences to introduce the K->A changes. Linking to an already linked residue
            is not allowed, and is caught as trying to change A->A.
        """
        for _, chB, resID in self.proper_links:

            resID = int(resID)

            if resID not in [6, 11, 27, 29, 33, 48, 63]:
                raise ValueError (f"Error, unsupported linkage {resID+1:}!")
            # Convert str resID to int and shift to 0-based
            resID = resID - 1

            pt = self.seq_ptr[chB]
            row = pt[0]
            col = pt[1] + resID
            if self.seq_tab[row][col] == 'K':
                seq = self.seq_tab[row]
                self.seq_tab[row] = seq[:col] + "A" + seq[col+1:]
            else:
                raise ValueError (f"Error, {chB}-{resID} was already mutated!")


    def create_sequence_list(self):

        # Get the valid row identifiers
        sequence_list = []
        for i,seq in enumerate(self.seq_tab):
            chainID = alpha[i]
            seq_rec = sequence_template(chainID, seq, use_template=self.use_template, use_msa=self.use_msa)
            sequence_list.append(seq_rec)

        # Create the linker instances. The C-ter G76-s are unique, and can be used as a linkerID.
        # To distinguish from chains, use a two-letter ID, starting with L.
        for linker in self.proper_links:
            print (f" Creating linker for {linker[0]}-{linker[1]}")
            # NOTE: this it the chainID after the lookup! Does not have to match linker[0] or linker[1]!
            chainID = alpha[self.seq_ptr[linker[0]][0]]
            sequence_list.append(linker_template("L"+chainID))

        # Extra sequences (non-ubiquitin) are prepended with an 'E'
        if self.extra_seq is not None:
            for i, seq in enumerate(args.extra_seq):
                seq_rec = sequence_template('E'+alpha[i],seq=seq)
                sequence_list.append(seq_rec)

        # Add sequences to the output json
        self.sequence_list = sequence_list
        self.out_json["sequences"] = sequence_list


    def create_bonds(self):
        # Create bonds section
        linker_list = []
        for cter_aa, mid_aa, link in self.proper_links:
            print (f" Linking {cter_aa}-G76 to {mid_aa}-K{link}")
            ca = alpha[self.seq_ptr[cter_aa][0]]
            ma = alpha[self.seq_ptr[mid_aa][0]]
            cterID = len(self.seq_tab[self.seq_ptr[cter_aa][0]])
            linker_list += bondedAtomPairs(ca,ma, cterID, int(link)+self.seq_ptr[mid_aa][1])

        self.out_json["bondedAtomPairs"] = linker_list



desc = """ This script adds ubiquitin chains with specific linkages to
    AF3 input. The linkers are denoted by two characters, starting
    with L: LA, LB, LC...

    Example input format:
       # B-G76 linked to A-K48...
       B A 48
       C B 48
       D C 48
       E A  6
       F -  0
"""

parser = argparse.ArgumentParser(
                    prog = 'Ubiquitin chain generator for AF3 predictions',
                    description = desc,
                    epilog = 'That\'s all folks!')
parser.add_argument('filename')           # positional argument
parser.add_argument('-n', '--name', type=str, help='name of the output json and the system', default='polyub')
parser.add_argument('--extra-seq', nargs="*", help='extra sequences')
parser.add_argument('--exclude-template', action='store_true', default=False)
parser.add_argument('--exclude-msa', action='store_true', default=False)
args = parser.parse_args()

ub = UbiquitinChain(args.filename,
                    name = args.name,
                    extra_seq = args.extra_seq,
                    use_template = not args.exclude_template,
                    use_msa = not args.exclude_msa)


if len(ub.m1_links) != 0:
    print ("Due to M1 linkages, the following chains have been relocated:")
    for key in ub.seq_ptr.keys():
        loc = ub.seq_ptr[key]
        print (key, "->", alpha[loc[0]], loc[1]+1)

