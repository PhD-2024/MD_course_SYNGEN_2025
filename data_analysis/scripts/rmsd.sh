#!/bin/bash

# Preprocessed trajectories (from preprocess.sh)
traj_list=("npt_ber_posres" "npt_ber" "npt" "nvt")

# Groups for RMSD 
groups=("DNA" "Backbone" "SideChain")

echo "=== Running RMSD analysis on centered trajectories ==="

for traj in "${traj_list[@]}"; do

    xtc="${traj}_center.xtc"
    tpr="${traj}_center.tpr"

    echo "Processing RMSD for: ${xtc}"

    for g in "${groups[@]}"; do

        out="rmsd_${traj}_g${g}.xvg"

        echo "  → RMSD group $g vs $g → $out"

        echo -e "${g}\n${g}" | gmx rms \
            -s "$tpr" \
            -f "$xtc" \
            -o "$out" \
            -tu ns

    done

    echo "Finished $traj"
    echo "--------------------------------------------------------"

done

echo "All RMSD analyses completed!"

