#!/bin/bash
source ./common.sh #to access set_common_vars() function
#The following function defines config, processes, version and datasets variables
set_common_vars "$1"
args=(
--config $config
--processes $processes #'dy_lep,vv,tt,st,wj,data'
--datasets $datasets
--version $version
--categories 'cat_mutau,mutau_signal_reg'
--cf.CalibrateEvents-workflow $workflow
--cf.SelectEvents-workflow $workflow
--cf.ReduceEvents-workflow $workflow
--cf.MergeReducedEvents-workflow local
--variables 'mutau_lep0_pt,mutau_lep0_eta,mutau_lep0_phi,mutau_lep0_ip_sig,'`
`'mutau_lep1_pt,mutau_lep1_eta,mutau_lep1_phi,mutau_lep1_mass,mutau_lep1_decayModePNet,mutau_lep1_decayMode,'`
`'mutau_mt,mutau_mvis,mutau_delta_r,mutau_pt,puppi_met_pt,puppi_met_phi' 
# 'etau_lep0_pt,etau_lep0_eta,etau_lep0_phi,etau_lep0_ip_sig,'`
# `'etau_lep1_pt,etau_lep1_eta,etau_lep1_phi,etau_lep1_mass,etau_lep1_decayModePNet,etau_lep1_decayMode,'`
# `'etau_mt,etau_mvis,etau_delta_r,puppi_met_pt,puppi_met_phi' 
--general-settings "cms-label=pw"
"${@:2}"
)
echo law run cf.PlotVariables1D "${args[@]}"
law run cf.PlotVariables1D "${args[@]}"



# etau channel 

