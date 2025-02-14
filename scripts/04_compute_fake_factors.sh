#!/bin/bash
source ./common_run3.sh #to access set_common_vars() function
#The following function defines config, processes, version and datasets variables
set_common_vars "$1"
args=(
        --config $config 
        --datasets $datasets
        --processes $processes
        --cf.CalibrateEvents-workflow $workflow
        --cf.CalibrateEvents-version $version
        
        --cf.SelectEvents-workflow $workflow
        --cf.SelectEvents-version $version

        --cf.ReduceEvents-workflow $workflow
        --cf.ReduceEvents-version $version
        
        --cf.MergeReducedEvents-workflow $workflow
        --cf.MergeReducedEvents-version $version

        --cf.ProvideReducedEvents-version $version

        --cf.ProduceColumns-producer 'ff_method'
         
        --cf.PrepareFakeFactorHistograms-version ff_method_fine_pt_binning
        --cf.ComputeFakeFactors-version ff_method_fine_pt_binning
       
        --cf.ComputeFakeFactors-categories 'cat_mutau_sr'
        "${@:2}"
    )
echo law run cf.ComputeFakeFactors "${args[@]}"
law run cf.ComputeFakeFactors "${args[@]}" 
