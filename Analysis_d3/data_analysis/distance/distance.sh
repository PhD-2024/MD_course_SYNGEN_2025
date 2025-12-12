gmx distance -s npt_center.tpr -f npt_center.xtc -select 'com of group "Protein" plus com of group "DNA"' -oall com_protein_dna.xvg
gmx distance -s npt.tpr -f npt.xtc -select 'com of group "Protein" plus com of group "DNA"' -oall com_protein_dna_before_proc.xvg

