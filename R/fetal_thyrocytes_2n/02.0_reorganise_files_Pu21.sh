#!/bin/bash

# Loop through all relevant files
for file in Data/published_scRNAseq/Pu_etal_2021/*.tsv.gz Data/published_scRNAseq/Pu_etal_2021/*.mtx.gz; do
    # Extract filename from full path
    filename=$(basename "$file")

    # Extract sample ID: second and third underscore-separated fields
    sample_id=$(echo "$filename" | cut -d'_' -f2,3)

    # Create destination directory (inside the same base path)
    dest_dir="Data/published_scRNAseq/Pu_etal_2021/$sample_id"
    mkdir -p "$dest_dir"

    # Determine new filename
    if [[ "$filename" == *barcodes.tsv.gz ]]; then
        newname="barcodes.tsv.gz"
    elif [[ "$filename" == *features.tsv.gz ]]; then
        newname="features.tsv.gz"
    elif [[ "$filename" == *matrix.mtx.gz ]]; then
        newname="matrix.mtx.gz"
    else
        echo "Unknown file type: $filename"
        continue
    fi

    # Move and rename the file
    mv "$file" "$dest_dir/$newname"
done
