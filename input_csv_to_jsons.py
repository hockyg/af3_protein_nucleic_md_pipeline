import os
import sys

def create_json_file(input_csv_row, output_json_file, mutid=0, ref=False):
    name = input_csv_row["uniprot_id"].lower()
    protein_seq = input_csv_row["peptide_sequence"]

    output_energy_file = output_json_file.replace(".input.json",".energy.txt")
    label = os.path.basename(output_json_file.replace(".input.json",""))

    if ref is True:
        nuc_strand_1 = input_csv_row["wild_nucleotide_strand1"]
        nuc_strand_2 = input_csv_row["wild_nucleotide_strand2"]
        dG = input_csv_row["wild_G0"]
        ddG = 0
    else:
        nuc_strand_1 = input_csv_row["mutated_nucleotide_strand1"]
        nuc_strand_2 = input_csv_row["mutated_nucleotide_strand2"]
        dG = input_csv_row["mutated_G0"]
        ddG = input_csv_row["ddG"]
    
    # write energy file
    with open(output_energy_file,'w') as fh:
        #add header
        fh.write(f"# name mutid dG ddG \n")
        fh.write(f"{name} {mutid} {dG} {ddG}\n")

    #get nuc name in lowercase
    nuc_name = input_csv_row["nucleotide_sequence_type"].lower()

    if(nuc_strand_2):
        print("2 nucleotide strands")
        num_nuc_chains = 2
        nuc_string = """{
                "%s": {
                "id": ["B"],
                "sequence": "%s"
                }
            },
            {
                "%s": {
                "id": ["C"],
                "sequence": "%s"
                }
            }"""%(nuc_name, nuc_strand_1.strip(), nuc_name, nuc_strand_2.strip())
    else:
        num_nuc_chains = 1
        nuc_string = """{%s: {
            "id": ["B"],
            "sequence": "%s"
            }}"""%(nuc_name, nuc_strand_1.strip()) 

    ion_chains = ",".join( numbers_to_chains(list(range(num_nuc_chains+2,num_nuc_chains+2+len(nuc_strand_1)))) )
    

    output_json_template = """{
    "name": "%s",
    "sequences": [
        {
            "protein": {
                "id": ["A"],
                "sequence": "%s"
            }
        },
        %s,
        {
            "ligand": {
                    "id": [%s],
                    "ccdCodes": ["MG"]
                }
        }
    ],
    "modelSeeds": [1],
    "dialect": "alphafold3",
    "version": 1
    }"""

    output_json = output_json_template%(label, protein_seq.strip(), nuc_string, ion_chains)

    print(f"# Writing sequence data for {nuc_name} + mut {mutid}, to {output_json_file}")
    with open(output_json_file,'w') as fh:
        fh.write(output_json)

def numbers_to_chains(numbers):
    def num_to_alpha(n):
        """Convert a number to its corresponding alphabetical sequence."""
        result = ""
        while n > 0:
            n -= 1  # Adjust for 0-based indexing
            result = chr(n % 26 + ord('A')) + result
            n //= 26
        return result

    return [f'"{num_to_alpha(num)}"' for num in numbers]


script=sys.argv[0]

if len(sys.argv)<2:
    print(f"Usage: python {script} input_csv")
    sys.exit()

input_file=sys.argv[1]
output_folder = os.path.join("output",os.path.basename(input_file).replace(".csv",""))
os.makedirs(output_folder,exist_ok=True)

#read csv file
import pandas as pd
df = pd.read_csv(input_file,keep_default_na=False)


#example
#              protein_name uniprot_id                                   peptide_sequence  ... wild_G0 mutated_G0   ddG
# 0  Cold shock protein CspB     P32081  MLEGKVKWFNSEKGFGFIEVEGQDDVFVHFSAIQGEGFKTLEEGQA...  ...   -8.55      -8.07  0.48
#headers: protein_name,uniprot_id,peptide_sequence,nucleotide_sequence_type,wild_nucleotide_strand1,wild_nucleotide_strand2,mutated_nucleotide_strand1,mutated_nucleotide_strand2,wild_G0,mutated_G0,ddG

name_count_dict = {}

json_files_to_run = []

for index, row in df.iterrows():
    name = row["uniprot_id"].lower()
    protein_output_dir = os.path.join(output_folder,name)
    os.makedirs(protein_output_dir,exist_ok=True)

    if name not in name_count_dict:
        name_count_dict[name] = 1
    else:
        name_count_dict[name] += 1

    mutid = name_count_dict[name]

    output_ref_json_file = os.path.join(protein_output_dir,f"{name}_ref.input.json")
    output_mut_json_file = os.path.join(protein_output_dir,f"{name}_mut{mutid}.input.json")

    if name_count_dict[name] == 1:
        json_files_to_run.append(output_ref_json_file)
    json_files_to_run.append(output_mut_json_file)

    if not os.path.exists(output_ref_json_file):
        create_json_file(row, output_ref_json_file, ref=True)
    if not os.path.exists(output_mut_json_file):
        create_json_file(row, output_mut_json_file, mutid=mutid, ref=False)

# print list of json files to a file in the output folder
with open(os.path.join(output_folder,"json_files_to_run.txt"),'w') as fh:
    fh.write("\n".join(json_files_to_run))
