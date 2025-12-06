echo -e "Protein\nDNA" | gmx hbond \
    -s npt_center.tpr \
    -f npt_center.xtc \
    -life hb_lifetime.xvg \
    -num hb_protein_dna.xvg 
