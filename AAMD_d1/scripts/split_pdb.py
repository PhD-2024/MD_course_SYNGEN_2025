''' this splits pdb into individual files '''
filename="files_afternoon_1/1J46.pdb"
with open(filename) as f:
    output_dict={}
    dump=[]
    part=0
    for line in f:
        if line.startswith("ATOM"):
            dump.append(line) 
        elif line.startswith("TER"):
            output_dict[part]=dump
            dump=[]
            part+=1
    for key in output_dict:
        with open(f"part_{key}.pdb","w") as out:
            out.writelines(output_dict[key])