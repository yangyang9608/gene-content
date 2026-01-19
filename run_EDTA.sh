#!/bin/bash

# Activate conda environment (if needed)
# conda activate /data/home/yangy/.conda/envs/EDTA

# Define list of species prefixes
species_list=("Adig" "Agra" "Agre" "Amad" "Aper" "Arub" "Asua" "Aza" "Bcei" "Cpen" "Och")

# Loop through each species (sequential execution)
for species in "${species_list[@]}"
do
    echo "========================================"
    echo "Processing species: $species"
    echo "Time: $(date)"
    echo "========================================"

    # Check if required input files exist
    if [[ ! -f "${species}.fa" ]]; then
        echo "Error: ${species}.fa file not found, skipping this species"
        continue
    fi

    if [[ ! -f "${species}.cds" ]]; then
        echo "Warning: ${species}.cds file not found, proceeding without CDS annotation"
        cds_param=""
    else
        cds_param="--cds ${species}.cds"
    fi
    
    # Run EDTA analysis (sequential execution, waiting for completion)
    echo "Starting EDTA analysis (will wait for completion)..."
    perl /home/yurm/software/EDTA/EDTA.pl \
        --species others \
        --genome ${species}.fa \
        ${cds_param} \
        --overwrite 0 \
        --sensitive 1 \
        --anno 1 \
        --force 1 \
        --threads 30 \
        2>&1 | tee ${species}_EDTA.log

    # Check if EDTA analysis succeeded
    if [ ${PIPESTATUS[0]} -eq 0 ]; then
        echo "EDTA analysis completed successfully for species $species"
    else
        echo "Error: EDTA analysis failed for species $species"
    fi

    echo "----------------------------------------"
done

echo "EDTA analysis completed for all species"
echo "Completion time: $(date)"
