import os
import sys

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

if len(sys.argv)<3:
    print(f"Usage: python {script} input_file output_file")
    sys.exit()

input_file=sys.argv[1]
output_json_file=sys.argv[2]
name=os.path.splitext( os.path.basename(input_file) )[0]

os.makedirs(os.path.dirname(output_json_file),exist_ok=True)

with open(input_file,'r') as fh:
    protein_seq, nuc_seq = fh.readlines()

if nuc_seq.find('U')>=0:
    nuc_name = "rna"
elif nuc_seq.find('T')>=0:
    nuc_name = "dna"
else:
    print("Error, could not determine nucleic acid type")
    sys.exit(1)
nuc_chains = ",".join( numbers_to_chains(list(range(3,3+len(nuc_seq)))) )

#output_json = template.format(protein_seq=protein_seq, nuc_seq=nuc_seq, nuc_chains=nuc_chains, name=name)
output_json_template = """{
  "name": "%s",
  "sequences": [
        {
          "protein": {
            "id": ["A"],
            "sequence": "%s"
          }
       },
       {
           "%s": {
                  "id": ["B"],
                  "sequence": "%s"
           }
       },
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

output_json = output_json_template%(name, protein_seq.strip(), nuc_name, nuc_seq.strip(), nuc_chains)

print(f"# Converting sequence data in {input_file} to {output_json_file}")
with open(output_json_file,'w') as fh:
    fh.write(output_json)

