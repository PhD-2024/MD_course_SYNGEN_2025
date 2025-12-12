#!/bin/bash

# List of EDR files to process
edr_list=("npt_ber.edr" "npt_ber_posres.edr" "npt.edr" "nvt.edr")

echo "=== Extracting Pressure, Density, Total Energy from EDR files ==="

for edr in "${edr_list[@]}"; do
    base="${edr%.edr}"
    echo "Processing: $edr"

    # --- Density ---
    echo "Density" | gmx energy -f "$edr" -o "${base}_density.xvg"

    # --- Total Energy or Potential (whichever exists) ---
    # First try Total Energy
    echo "Total Energy" | gmx energy -f "$edr" -o "${base}_total_energy.xvg"

    # If Total Energy wasn't found, use Potential instead
    if [ ! -s "${base}_total_energy.xvg" ]; then
        echo "Potential" | gmx energy -f "$edr" -o "${base}_total_energy.xvg"
    fi

    echo " ${base}_density.xvg, ${base}_total_energy.xvg"
    echo "--------------------------------------------------------------"

done

echo "Energy extraction complete!"

