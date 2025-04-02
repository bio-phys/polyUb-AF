import json
import numpy as np
import argparse
from string import ascii_uppercase


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
        self.ub_ids = ub_ids[~mask]
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
        # Create an Ub chain table and Ub chain pointers
        seq_tab = []
        seq_ptr = {}

        for i, ub_id in enumerate(self.ub_ids):
            seq_tab.append(UbiquitinChain.ub_seq)
            seq_ptr[ub_id] = [i,0]

        # Introduce all the M1 linkages
        for chA, chB, resID in self.m1_links:
            ptA = seq_ptr[chA]
            ptB = seq_ptr[chB]
            # Elongate the chain
            seq_tab[ptB[0]] += self.ub_seq
            # Shift all indices in this row
            for key in seq_ptr.keys():
                if seq_ptr[key][0] == ptB[0]:
                    seq_ptr[key][1] += len(self.ub_seq)
            # Update the first pointer
            ptA[0] = ptB[0]

        self.seq_tab = seq_tab
        self.seq_ptr = seq_ptr
        # Only these rows of seq_tab are valid
        self.valid_rows = np.sort(np.unique([elem[0] for elem in seq_ptr.values()]))


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
        for i,r in enumerate(self.valid_rows):
            seq_rec = sequence_template(ascii_uppercase[i], self.seq_tab[r], use_template=self.use_template, use_msa=self.use_msa)
            sequence_list.append(seq_rec)

        # Create the linker instances. The C-ter G76-s are unique, and can be used as a linkerID.
        # To distinguish from chains, use a two-letter ID, starting with L.
        for linker in self.proper_links:
            print (f" Creating linker for {linker[0]}-{linker[1]}")
            sequence_list.append(linker_template("L"+linker[0]))

        # Extra sequences (non-ubiquitin) are prepended with an 'E'
        if self.extra_seq is not None:
            for i, seq in enumerate(args.extra_seq):
                seq_rec = sequence_template('E'+ascii_uppercase[i],seq=seq)
                sequence_list.append(seq_rec)

        # Add sequences to the output json
        self.sequence_list = sequence_list
        self.out_json["sequences"] = sequence_list


    def create_bonds(self):
        # Create bonds section
        linker_list = []
        for cter_aa, mid_aa, link in self.proper_links:
            print (f" Linking {cter_aa}-G76 to {mid_aa}-K{link}")
            ca = ascii_uppercase[self.seq_ptr[cter_aa][0]]
            ma = ascii_uppercase[self.seq_ptr[mid_aa][0]]
            cterID = len(self.seq_tab[self.seq_ptr[cter_aa][0]])
            linker_list += bondedAtomPairs(ca,ma, cterID, int(link)+self.seq_ptr[ma][1])

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
parser.add_argument('-n', '--name', type=str, help='name of the output json and the system')
parser.add_argument('--extra-seq', nargs="*", help='extra sequences')
parser.add_argument('--exclude-template', action='store_true', default=False)
parser.add_argument('--exclude-msa', action='store_true', default=False)
args = parser.parse_args()

# print (args)

ub = UbiquitinChain(args.filename,
                    name = args.name,
                    extra_seq = args.extra_seq,
                    use_template = not args.exclude_template,
                    use_msa = not args.exclude_msa)
