#!/bin/bash

set_common_vars() {

version="fixed_also_mutau_trigmatch"
case $1 in
    "run3ts_lim")
        config="run3_2022_preEE_tau_spinner_limited"
        datasets='h_ggf_htt_filtered' #,h_ggf_htt_unfiltered,zh_htt_unfiltered'
        processes='h_ggf_htt,data_mu' #,zh_htt'
        workflow='local'
    ;;
    "run3ts")
        config="run3_2022_preEE_tau_spinner"
        datasets='h_ggf_htt_filtered' #,h_ggf_htt_unfiltered,zh_htt_unfiltered'
        processes='h_ggf_htt' #,zh_htt'
        workflow='htcondor'
    ;;
    "run3_lim")
        config="run3_2022_preEE_limited"
        datasets='dy_incl' 
        processes='dy_z2mumu,dy_lep' 
        workflow='local'
    ;;
    "run3_data")
        config="run3_2022_preEE"
        datasets='data_mu_C,data_mu_D,data_tau_C,data_tau_D,data_e_C,data_e_D'
        processes='data' 
        workflow='htcondor'
    ;;
    "run3_data_lim")
        config="run3_2022_preEE_limited"
        datasets='data_e_C'
        processes='data' 
        workflow='local'
    ;;
    "run3preEE")
        config="run3_2022_preEE"
        data='data_mu_C,data_mu_D,data_tau_C,data_tau_D,data_e_C,data_e_D,'
        bkg_ewk='wj_incl,ww,wz,zz,dy_incl,'
        bkg_top='st_twchannel_t_fh,st_twchannel_t_sl,st_twchannel_t_dl,st_twchannel_tbar_sl,st_twchannel_tbar_dl,st_tchannel_tbar,st_tchannel_t,'
        bkg_ttbar='tt_sl,tt_dl,tt_fh'
        datasets="$data$bkg_ewk$bkg_top$bkg_ttbar"
        processes='dy_z2tautau,dy_z2mumu,dy_z2ee,vv,tt,st,wj,data'
        workflow='htcondor'
    ;;
    *)
    echo "Unknown run argument!"
    exit
    ;;
esac
}
