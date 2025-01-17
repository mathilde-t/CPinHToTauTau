#!/bin/bash

set_common_vars() {

version="dev"
    
categories_mutau="mutau_signal_reg,mutau_signal_reg_no_mt,mutau_signal_reg_endcap_tau,mutau_signal_reg_barrel_tau,mutau_signal_reg_endcap_tau_no_mt,mutau_signal_reg_barrel_tau_no_mt" 
variables_mutau='mutau_lep0_pt,mutau_lep0_eta,mutau_lep0_phi,mutau_lep0_ip_sig,mutau_lep1_pt,mutau_lep1_eta,mutau_lep1_phi,mutau_lep1_mass,mutau_lep1_decayModePNet,mutau_lep1_decayMode,mutau_mt,mutau_mvis,mutau_delta_r,mutau_pt,puppi_met_pt,puppi_met_phi'

categories_etau="etau_signal_reg,etau_signal_reg_no_mt,etau_signal_reg_endcap_tau,etau_signal_reg_barrel_tau,etau_signal_reg_endcap_tau_no_mt,etau_signal_reg_barrel_tau_no_mt" 
variables_etau='jet_1_pt,etau_lep0_pt,etau_lep0_eta,etau_lep0_phi,etau_lep0_ip_sig,etau_lep1_pt,etau_lep1_eta,etau_lep1_phi,etau_lep1_mass,etau_lep1_decayModePNet,etau_lep1_decayMode,etau_mt,etau_mvis,etau_delta_r,etau_pt,puppi_met_pt,puppi_met_phi'

data_e_2022preEE='data_e_C,data_e_D,'
data_mu_2022preEE='data_singlemu_C,data_mu_C,data_mu_D,'

data_e_2022postEE='data_e_E,data_e_F,data_e_G,'
data_mu_2022postEE='data_mu_E,data_mu_F,data_mu_G,'

bkg_ewk='wj_incl_madgraph,ww,wz,zz,dy_lep_madgraph,'
bkg_top='st_twchannel_t_sl,st_twchannel_t_dl,st_twchannel_tbar_sl,st_twchannel_tbar_dl,st_tchannel_tbar,st_tchannel_t,st_schannel_t_lep,st_schannel_tbar_lep,'
bkg_ttbar='tt_sl,tt_dl,tt_fh'

data_e_2023preBPix='data_e_Cv123,data_e_Cv4'
data_e_2023postBPix='data_e_D'
data_mu_2023preBPix='data_m_Cv123,data_m_Cv4'
data_mu_2023postBPix='data_m_D'


case $1 in
##############################
####### 2022preEE ############
##############################
    "run3_2022preEE_etau_lim")
        config="run3_2022_preEE_etau_limited"	
        datasets='ww'
        processes='dy_z2mumu,dy_z2ee,dy_lep,vv'
	categories=$categories_etau
	variables=$variables_etau
        workflow='htcondor'
    ;;
    "run3_2022preEE_etau")
        config="run3_2022_preEE_etau"
        data=$data_e_2022preEE
        bkg_ewk=$bkg_ewk
        bkg_top=$bkg_top
        bkg_ttbar=$bkg_ttbar
        datasets="$data$bkg_ewk$bkg_top$bkg_ttbar"
        processes='dy_z2tautau,dy_z2mumu,dy_z2ee,vv,tt,st,wj,data'
	categories=$categories_etau
	variables=$variables_etau
        workflow='htcondor'
    ;;
    "run3_2022preEE_mutau")
        config="run3_2022_preEE_mutau"
        data=$data_mu_2022preEE
        bkg_ewk=$bkg_ewk
        bkg_top=$bkg_top
        bkg_ttbar=$bkg_ttbar
        datasets="$data$bkg_ewk$bkg_top$bkg_ttbar"
        processes='dy_z2tautau,dy_z2mumu,dy_z2ee,vv,tt,st,wj,data'
	categories=$categories_mutau
	variables=$variables_mutau
        workflow='htcondor'
     ;;
###############################
######## 2022postEE ###########
###############################    
    "run3_2022postEE_etau_lim")
        config="run3_2022_postEE_etau_limited"	
        datasets='data_e_E' 
        processes='data' 
	categories=$categories_etau
	variables=$variables_etau
        workflow='local'
    ;;
    "run3_2022postEE_mutau_lim")
        config="run3_2022_postEE_mutau_limited"	
        datasets='dy_lep_madgraph' 
        processes='dy_z2mumu,dy_z2ee,dy_lep' 
	categories=$categories_mutau
	variables=$variables_mutau
        workflow='local'
    ;;
    "run3_2022postEE_etau")
        config="run3_2022_postEE_etau"
        data=$data_e_2022postEE
        bkg_ewk=$bkg_ewk
        bkg_top=$bkg_top
        bkg_ttbar=$bkg_ttbar
        datasets=$data$bkg_ewk$bkg_top$bkg_ttbar
        processes='dy_z2tautau,dy_z2mumu,dy_z2ee,vv,tt,st,wj,data'
	    categories=$categories_etau
	    variables=$variables_etau
	    workflow='htcondor'
    ;;
    "run3_2022postEE_mutau")
        config="run3_2022_postEE_mutau"
        data=$data_mu_2022postEE
        bkg_ewk=$bkg_ewk
        bkg_top=$bkg_top
        bkg_ttbar=$bkg_ttbar
	datasets="$data$bkg_ewk$bkg_top$bkg_ttbar"
        processes='dy_z2tautau,dy_z2mumu,dy_z2ee,vv,tt,st,wj,data'
	categories=$categories_mutau
	variables=$variables_mutau
	workflow='htcondor'
    ;;
    *)
    echo "Unknown run argument!"
    exit
    ;;
##############################
####### 2023preBPix ############
##############################

    "run3_2023preBPix_etau")
        config="run3_2023_preBPix_etau"
        data=$data_e_2023preBPix
        bkg_ewk=$bkg_ewk
        bkg_top=$bkg_top
        bkg_ttbar=$bkg_ttbar
        datasets="$data$bkg_ewk$bkg_top$bkg_ttbar"
        processes='dy_z2tautau,dy_z2mumu,dy_z2ee,vv,tt,st,wj,data'
    categories=$categories_etau
    variables=$variables_etau
        workflow='htcondor'
    ;;
    "run3_2023preBPix_mutau")
        config="run3_2023_preBPix_mutau"
        data=$data_mu_2023preBPix
        bkg_ewk=$bkg_ewk
        bkg_top=$bkg_top
        bkg_ttbar=$bkg_ttbar
        datasets="$data$bkg_ewk$bkg_top$bkg_ttbar"
        processes='dy_z2tautau,dy_z2mumu,dy_z2ee,vv,tt,st,wj,data'
    categories=$categories_mutau
    variables=$variables_mutau
        workflow='htcondor'
     ;;


##############################
####### 2023postBPix ############
##############################

    "run3_2023postBPix_etau")
        config="run3_2023_postBPix_etau"
        data=$data_e_2023postBPix
        bkg_ewk=$bkg_ewk
        bkg_top=$bkg_top
        bkg_ttbar=$bkg_ttbar
        datasets="$data$bkg_ewk$bkg_top$bkg_ttbar"
        processes='dy_z2tautau,dy_z2mumu,dy_z2ee,vv,tt,st,wj,data'
    categories=$categories_etau
    variables=$variables_etau
        workflow='htcondor'
    ;;
    "run3_2023postBPix_mutau")
        config="run3_2023_postBPix_mutau"
        data=$data_mu_2023postBPix
        bkg_ewk=$bkg_ewk
        bkg_top=$bkg_top
        bkg_ttbar=$bkg_ttbar
        datasets="$data$bkg_ewk$bkg_top$bkg_ttbar"
        processes='dy_z2tautau,dy_z2mumu,dy_z2ee,vv,tt,st,wj,data'
    categories=$categories_mutau
    variables=$variables_mutau
        workflow='htcondor'
     ;;


esac
}
