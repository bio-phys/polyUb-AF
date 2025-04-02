import json
import sys
import numpy as np

""" This script add ubiquitin chains with specific linkages to
    AF3 input. The linkers are denoted by two characters, starting
    with L: LA, LB, LC...
"""

def sequence_template(seq_id, use_template=True, use_msa=True):
    """ Print a single wild-type ubiquitin chain with a specific chainID
    """

    wt_ubiquitin = "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"

    seq = {
      "protein": {
        "id": seq_id,
        "sequence": wt_ubiquitin,
      }
    }

    # Set explicitly to empty string if not used
    if not use_template:
        seq["protein"]["templates"] = ""

    if not use_msa:
        seq["protein"]["unpairedMsa"] = ""
        seq["protein"]["pairedMsa"] = ""

    return seq


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

def bondedAtomPairs(start, end, link, linker_A="C1", linker_B="C3"):
    """
    """
    return [
      [[end  , int(link),  "CB"], ["L"+start, 1, linker_A]],
      [[start,        76, "OXT"], ["L"+start, 1, linker_B]]
    ]

# The json dictionary to populate.
out_json = {
    "name": "polyub",
    "modelSeeds": [1,],
    "dialect": "alphafold3",
    "version": 2,
    "sequences": [],
}

# Read the file with Ub linkages
# Example format:
#    # B-G76 linked to A-K48...
#    B A 48
#    C B 48
#    D C 48
#    E A  6
links = np.loadtxt(sys.argv[1], dtype=str)

# Get the unique Ub chains
ub_ids = np.unique(links[:,:2].flatten())
print (f"Found {len(ub_ids)} unique ubiquitin chains in {len(links)} linkages!")

# Create the Ub instances
sequence_list = []
for ub_id in ub_ids:
    sequence_list.append(sequence_template(ub_id))

# Modify the Ub sequences to introduce the K->A changes. Linking to an already linked residue
# is not allowed, and is caught as trying to change A->A.
for chainID, resID in links[:,1:]:

    print (f" Mutating {chainID}-{resID}")

    # Convert str resID to int and shift to 0-based
    resID = int(resID) - 1
    for sequence in sequence_list:
        if sequence['protein']['id'] == chainID:

            # Short-hand for the sequence string
            seq = sequence['protein']['sequence']

            # Check if not already an Alanine
            if seq[resID] == "A":
                raise ValueError (f"Error, {chainID}-{resID} was already mutated!")
            else:
                # Mutate to A
                sequence['protein']['sequence'] = seq[:resID] + "A" + seq[resID+1:]

# Create the linker instances. The C-ter G76-s are unique, and can be used as a linkerID.
# To distinguish from chains, use a two-letter ID, starting with L.
for linker in links[:,0]:
    sequence_list.append(linker_template("L"+linker))

# Add sequences to the output json
out_json["sequences"] = sequence_list

# Create bonds section
linker_list = []
for cter_aa, mid_aa ,link in links:
    print (f" Linking {cter_aa}-G76 to {mid_aa}-K{link}")
    linker_list += bondedAtomPairs(cter_aa,mid_aa,link)

out_json["bondedAtomPairs"] = linker_list


with open('data.json', 'w', encoding='utf-8') as f:
    json.dump(out_json, f, indent=4, ensure_ascii=False)

