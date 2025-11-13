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


echo "1" | gmx pdb2gmx -f part_2_copy.pdb -ff amber99bsc1 -ignh -o protein.gro
