import argparse
parser=argparse.ArgumentParser(description="Convert files between formats.")
parser.add_argument("--file", required=True )
args=parser.parse_args()
with open(args.file,"r") as inf:
	printing=False
	for line in inf:
		if "[ moleculetype ]" in line:
			printing=True
		if "; Include water topology" in line:
			printing=False

		if printing and (len(line.split(sep=";")[0])>= 1):
			if "helix_part/posre" in line:
				line=line.replace("helix_part/posre","posre")
			elif "protein_part/posre" in line:
				line=line.replace("protein_part/posre","posre")
			print(line.split(sep=";")[0])

