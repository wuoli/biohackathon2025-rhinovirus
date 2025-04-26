#!/bin/bash

#Note: To activate EDirect for this terminal session, please execute the following:
#export PATH=${HOME}/edirect:${PATH}

# Ensure output file is empty before starting
> westnile_metadata.tsv

for id in $(cat westnile_accessions_cleaned.txt); do
    echo "Processing $id"

    # Fetch the GenBank entry
    entry=$(esearch -db nuccore -query "$id" | efetch -format gb)

    # Extract geo_loc_name and collection_date using awk
    geo_loc=$(echo "$entry" | awk -F'"' '/\/geo_loc_name=/ {print $2; exit}')
    collection_date=$(echo "$entry" | awk -F'"' '/\/collection_date=/ {print $2; exit}')

    # Handle missing fields
    geo_loc=${geo_loc:-NA}
    collection_date=${collection_date:-NA}

    # Append to metadata.tsv
    echo -e "$id\t$geo_loc\t$collection_date" >> westnile_metadata.tsv
done
