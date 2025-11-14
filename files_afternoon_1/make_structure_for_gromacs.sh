#First we split our pdb
python3 split_pdb.py
#last part was the protein 
#if we convet it directly to a topology, we get issues because naming was different for the different expectations
#so we rename it

cp part_2.pdb part_2_copy.pdb
line_to_replace=$(grep "HB3 MET A   1" part_2.pdb)
replace_by=$(sed "s|HB3|HB2|g" <<< ${line_to_replace})

echo "names do not fit together"
echo replacing $line_to_replace
echo by $replace_by

#also another case

line_to_replace=$(grep "HG3 MET A   1" part_2.pdb)
replace_by=$(sed "s|HG3|HG2|g" <<< ${line_to_replace})


sed -i "s|${line_to_replace}|${replace_by}|g" part_2_copy.pdb

mkdir -p protein_part
mkdir -p helix_part

echo "1" | gmx pdb2gmx -f part_0.pdb -ff amber99bsc1 -ignh -o  helix_part/helix_pt1.gro -p helix_part/helix_topol_pt1.top -i helix_part/posre_h_pt1.itp
echo "1" | gmx pdb2gmx -f part_1.pdb -ff amber99bsc1 -ignh -o  helix_part/helix_pt2.gro -p helix_part/helix_topol_pt2.top -i helix_part/posre_h_pt2.itp
echo "1" | gmx pdb2gmx -f part_2_copy.pdb -ff amber99bsc1 -ignh -o protein_part/protein.gro -p protein_part/protein_topol.top -i protein_part/posre.itp



#we want 3 GEOMETRY files for the MD-simulation 

# both helix strands 
# protein part
# the whole system


#so far we have all 3 Molecules (strand 1 strand2 and the protein part separately)
#putting them together is simple:
#we just need to make a single .gro file.
#a gro file has a header- here we do not care what it says.
#then there is the number of atoms- that we do need to modify
#
#then there is a coordinate block-
#for this take the individual ones and copy the coordinate lines
#at the end we put the boxdefintion again

#number of lines in the helix_pt1.gro
nh1=$(cat helix_part/helix_pt1.gro | wc -l)
nh2=$(cat helix_part/helix_pt2.gro | wc -l)
np=$(cat protein_part/protein.gro | wc -l)
num_atoms=$(( $nh1 + $nh2 +$np - 9)) #3 extra line per file
num_atoms_H=$(( $nh1 + $nh2 -6 )) # same here

echo "full system" > full_system.gro
echo "helix system" > helix_system.gro
echo $num_atoms >> full_system.gro
echo $num_atoms_H >> helix_system.gro

head -n $(($nh1 -1)) helix_part/helix_pt1.gro | tail -n $(( $nh1 -3 )) >> helix_system.gro
head -n $(($nh1 -1)) helix_part/helix_pt1.gro | tail -n $(( $nh1 -3 )) >> full_system.gro
#here we also want the boxsize
tail -n $(( $nh2 -2 )) helix_part/helix_pt2.gro >> helix_system.gro
head -n $(($nh2 -1)) helix_part/helix_pt2.gro | tail -n $(( $nh2 -3 )) >> full_system.gro
#here we also want the boxsize
tail -n $(( $np -2 )) protein_part/protein.gro >> full_system.gro

#check  how large the structure is and modify the boxsize, so that the largest distance is AT LEAST half the box dimension
# why is that necessary?
topfiles="full_system.gro helix_system.gro"

for i in $topfiles
do
	what=$(tail -n1  $i)
	sed -i "s|$what|11.0 11.0 11.0|g" $i 
done

#we also want a topology file that can do this
#essentially a top file has the atomtype definitions at the beginning and then uses all of those for the moleculefiles
#we can use the  generated topologies and make a shared topfile for all

cat "./amber99bsc1.ff/forcefield.itp" > system.top
echo '#include "
