#!/bin/bash

# Trajectories to preprocess
traj_list=("npt") # for all files ("npt_ber_posres" "npt_ber" "npt" "nvt")

echo "=== Preprocessing trajectories (nojump → center + pbc mol, KEEPING WATER) ==="

for traj in "${traj_list[@]}"; do

    xtc="${traj}.xtc"
    tpr="${traj}.tpr"
    gro="${traj}.gro"
    echo "Processing $xtc"

    # Step 1: Remove PBC jumps (reconstruct molecules)
    
    echo -e "non-Water" | gmx trjconv -s "$tpr" -f "$xtc" -o "${traj}_nojump.xtc" -pbc nojump 
    echo -e "non-Water" |  gmx trjconv -s "$tpr" -f "$xtc" -o "${traj}_nojump.gro" -pbc nojump -b 0 -e 0
    
    gmx grompp -f npt.mdp -c "${traj}_nojump.gro" -o "${traj}_nojump.tpr" -p full_system_ions_dry.top -maxwarn 1   
    
    echo -e "System\nSystem" | gmx trjconv -s "${traj}_nojump.tpr" -f "${traj}_nojump.xtc" -o "${traj}_cluster.xtc" -pbc cluster
    echo -e "DNA\nSystem" | gmx trjconv -s "${traj}_nojump.tpr" -f "${traj}_cluster.xtc" -o "${traj}_center.xtc" -center
    
    echo -e "System\nSystem" | gmx trjconv -s "${traj}_nojump.tpr" -f "${traj}_center.xtc" -o "${traj}_center.gro" -b 0 -e 0
    
    gmx grompp -f npt.mdp -c "${traj}_center.gro" -o "${traj}_center.tpr" -p full_system_ions_dry.top -maxwarn 1
    echo "Clean (PBC-corrected) trajectory saved → ${traj}_nojump.xtc"
    echo "-------------------------------------------------------------"

done

echo "Preprocessing completed."

