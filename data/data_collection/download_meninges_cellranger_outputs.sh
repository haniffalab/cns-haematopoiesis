#!/bin/bash

function execute_iget_command() {
    local file_list="$1"
    local destination="$2"

    while IFS= read -r filepath; do
        iget -r -KVf --progress "$filepath" "$destination"
    done < "$file_list"
}

# Usage example:
execute_iget_command "/lustre/scratch126/cellgen/team298/sm54/Data_Integration/Meninges/src/cellranger_files_iRODS_paths_meninges.txt" "/lustre/scratch126/cellgen/team298/sm54/Data_Integration/Meninges/data/cellranger_outputs_all_samples/"