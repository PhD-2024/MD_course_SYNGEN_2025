echo -e 'Protein\n' | gmx rmsf -s npt_center.tpr -f npt_center.xtc -o rmsf_protein.xvg -res
echo -e 'DNA\n' | gmx rmsf -s npt_center.gro -f npt_center.xtc -o rmsf_DNA.xvg -res


