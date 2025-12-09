echo -e 'Protein\nDNA' | gmx mindist -s npt_center.tpr -f npt_center.xtc -group -od mindist_protein_dna.xvg

