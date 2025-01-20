#!/bin/bash
source ./common_run3.sh #to access set_common_vars() function
#The following function defines config, processes, version and datasets variables
set_common_vars "$1"
args=(
        --config $config
        --datasets dy_lep_madgraph,wj_incl_madgraph
        --processes dy_z2tautau,dy_z2mumu,dy_z2ee,wj
        --categories 'mutau_signal_reg'
        --version dev_deepTau-MTT_nob_veto

        --cf.CalibrateEvents-workflow $workflow
        --cf.CalibrateEvents-version $version
        
        --cf.SelectEvents-workflow $workflow
        --cf.SelectEvents-version $version
        
        --cf.MergeReducedEvents-workflow $workflow
        --cf.MergeReducedEvents-version $version

        --cf.ProvideReducedEvents-version $version
        --cf.CreateFakeFactorHistograms-weight-producer ff_weight_producer
        "${@:2}"
        
    )
echo law run cf.MergeFakeFactors "${args[@]}"
law run cf.MergeFakeFactors "${args[@]}" 