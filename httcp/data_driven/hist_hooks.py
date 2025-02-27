"""
Histogram hooks.
"""

from collections import defaultdict

import law
import order as od

import scinum as sn

from columnflow.util import maybe_import, DotDict

np = maybe_import("numpy")
hist = maybe_import("hist")


logger = law.logger.get_logger(__name__)


def add_hist_hooks(config: od.Config) -> None:
    """
    Add histogram hooks to a configuration.
    """
    # helper to convert a histogram to a number object containing bin values and uncertainties
    # from variances stored in an array of values
    def hist_to_num(h: hist.Hist, unc_name=str(sn.DEFAULT)) -> sn.Number:
        return sn.Number(h.values(), {unc_name: h.variances()**0.5})
    
    
    def qcd_estimation(task, hists, category_inst):
        
        def rel_err(h_arr=[], err_arr=[]):
            if len(h_arr):  sum_var = np.zeros_like(h_arr[0].values())
            else: sum_var = np.zeros_like(err_arr[0])
            for x in h_arr: sum_var += x.variances()/np.maximum(x.values()**2, 1)
            for the_arr in err_arr: sum_var += err_arr
            return np.sqrt(sum_var)
        
        def get_hists_from_reg(config: od.Config, hists: dict, region: str)-> hist.Hist :
            hists_ = hists[region]
            cat_id = config.get_category(region).id
            data_hists = []
            mc_hists = []
            for (proc, h) in hists_.items():
                if proc.is_data:
                    data_hists.append(h)
                elif proc.is_mc:
                    mc_hists.append(h)
            
            mc_hist = sum(mc_hists[1:], mc_hists[0].copy())
            data_hist = sum(data_hists[1:], data_hists[0].copy())
            
            return data_hist, mc_hist
        
        sr = category_inst
        data_num, mc_num = get_hists_from_reg(config, hists,sr.aux['abcd_regs']['dr_num'])
        data_den, mc_den = get_hists_from_reg(config, hists, sr.aux['abcd_regs']['dr_den']) 
        data_ar, mc_ar = get_hists_from_reg(config, hists,sr.aux['abcd_regs']['ar'])
        num = np.sum(data_num.values()) - np.sum(mc_num.values())
        den = np.sum(data_den.values()) - np.sum(mc_den.values()) 
        
        tf = num/np.maximum(den,0.001)
        tf_err2 = ((np.sum(data_num.variances()) + np.sum(mc_num.variances()))/den**2 + 
                  tf**2/den**2 *(np.sum(data_den.variances()) + np.sum(mc_den.variances())))
        
        from cmsdb.processes.qcd import qcd
        hists_sr = hists[sr.name].copy()
        h_donor_name = list(hists_sr.keys())[0]
        hists_sr[qcd] = hists_sr[h_donor_name].copy().reset()
        hists_sr[qcd].view().value = np.maximum(data_ar.values() -  mc_ar.values(), 0.) * tf
        #hists_sr[qcd].view().variance = data_ar.values()**2 * tf_err2 + data_ar.variances() * (tf**2)
        return hists_sr
    
    def ff_method(task, hists, category_inst):
        if not hists:
            return hists
        # extract all unique category ids and verify that the axis order is exactly
        # "category -> shift -> variable" which is needed to insert values at the end
        def rel_err(h_arr=[], err_arr=[]):
            if len(h_arr):  sum_var = np.zeros_like(h_arr[0].values())
            else: sum_var = np.zeros_like(err_arr[0])
            for x in h_arr: sum_var += x.variances()/np.maximum(x.values()**2, 1)
            for the_arr in err_arr: sum_var += err_arr
            return np.sqrt(sum_var)

        sig_reg = category_inst
        
        def get_ar_data_hist(config: od.Config, hists: dict, region: str)-> hist.Hist :
            h_reg = hists[region]
            data_hists = [h for p, h in h_reg.items() if p.is_data]
            data_hist = sum(data_hists[1:], data_hists[0].copy())
            cat_id = config.get_category(region).id
            return data_hist
        
        def calc_wj_yields(config: od.Config, hists: dict, region: str)-> hist.Hist :
            h_reg = hists[region]
            data_hists = [h for p, h in h_reg.items() if p.is_data]
            data_hist = sum(data_hists[1:], data_hists[0].copy())
            cat_id = config.get_category(region).id
            data_ar = data_hist
            
            wj_hists = [h for p, h in h_reg.items() if 'wj' in p.name]
            wj_hist = sum(wj_hists[1:], wj_hists[0].copy())
            wj_ar = wj_hist
            
            r = wj_ar.values() / np.maximum(data_ar.values(), 1)
            r_err =  np.where((data_ar.values() > 0),
                              rel_err(h_arr=[wj_ar,data_ar]) * r,
                              np.ones_like(r)) 
            return r[0], r_err[0]
        
        hist_qcd = get_ar_data_hist(config,hists, sig_reg.replace('sr', 'ar_qcd'))
        hist_wj  = get_ar_data_hist(config,hists, sig_reg.replace('sr', 'ar_wj'))
        wj_yield, yield_err = calc_wj_yields(config,hists, sig_reg.replace('sr', 'ar_yields'))
      
        fakes_val = hist_wj.values() * wj_yield + hist_qcd.values() * (1. - wj_yield)
        
        #Carefull this error estimation is not quite correct because yields errors and hist errors are corellated
        fakes_err = rel_err([hist_wj,hist_qcd],err_arr=[yield_err]) * fakes_val
        
        hists_sr = hists[sig_reg].copy()
        
        #Remove wj histogram from the signal region set
        from cmsdb.processes.qcd import jet_fakes
        
        wj_proc = [p for p in hists_sr.keys() if 'wj' in p.name]
        hists_sr[jet_fakes] = hists_sr[wj_proc[0]].copy().reset()
        cat_id = config.get_category(sig_reg).id
        hists_sr[jet_fakes].view().value = fakes_val
        hists_sr[jet_fakes].view().variance = fakes_err
        
        del hists_sr[wj_proc[0]]   
        return hists_sr
    

    config.x.hist_hooks = {
        "good_old_abcd"  : qcd_estimation,
        "ff_method" : ff_method,
    }
