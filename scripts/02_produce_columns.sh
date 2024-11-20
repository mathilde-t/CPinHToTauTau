#!/bin/bash
source ./common.sh #to access set_common_vars() function
#The following function defines config, processes, version and datasets variables
set_common_vars "$1"
args=(
        --config $config
        --datasets $datasets
        --version $version
        "${@:2}"
    )
    
echo law run cf.ProduceColumnsWrapper "${args[@]}"
law run cf.ProduceColumnsWrapper "${args[@]}"




