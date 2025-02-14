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

    # # helper to integrate values stored in an array based number object
    # def integrate_num(num: sn.Number, axis=None) -> sn.Number:
    #     return sn.Number(
    #         nominal=num.nominal.sum(axis=axis),
    #         uncertainties={
    #             unc_name: (
    #                 (unc_values_up**2).sum(axis=axis)**0.5,
    #                 (unc_values_down**2).sum(axis=axis)**0.5,
    #             )
    #             for unc_name, (unc_values_up, unc_values_down) in num.uncertainties.items()
    #         },
    #     )

    def qcd_estimation(task, hists):
        if not hists:
            return hists
        from cmsdb.processes.qcd import qcd
        # extract all unique category ids and verify that the axis order is exactly
        # "category -> shift -> variable" which is needed to insert values at the end
        CAT_AXIS, SHIFT_AXIS, VAR_AXIS = range(3)
        category_ids = set()
        for proc, h in hists.items():
            # validate axes
            assert len(h.axes) == 3
            assert h.axes[CAT_AXIS].name == "category"
            assert h.axes[SHIFT_AXIS].name == "shift"
            # get the category axis
            cat_ax = h.axes["category"]
            for cat_index in range(cat_ax.size):
                category_ids.add(cat_ax.value(cat_index))

        # create qcd groups
        signal_categories : dict[str, dict[str, od.Category]] = defaultdict(DotDict)
        #Filling the dictionay with signal histograms
        signal_category = config.get_category(task.categories[0])

        if signal_category.aux:

            # sum up mc and data histograms, stop early when empty
            mc_hists = [h for p, h in hists.items() if p.is_mc and not p.has_tag("signal")]
            data_hists = [h for p, h in hists.items() if p.is_data]
            mc_hist = sum(mc_hists[1:], mc_hists[0].copy())
            data_hist = sum(data_hists[1:], data_hists[0].copy())

            # start by copying the mc hist and reset it, then fill it at specific category slices
            hists[qcd] = qcd_hist = mc_hist.copy().reset()
            
            def get_hist (h, category): 
                return h[{"category": hist.loc(category.id)}]
            
            
            def find_by_id(hist, cat_id):
                cat_axis = hist.axes["category"]
                for  the_idx, the_id in enumerate(cat_axis): 
                    if the_id == cat_id: return the_idx
                return -1

            control_cat = config.get_category(signal_category.aux['control_reg'][0])
            dr_num = config.get_category(signal_category.aux['control_reg'][0].replace('ar','dr_num'))
            dr_den = config.get_category(signal_category.aux['control_reg'][0].replace('ar','dr_den'))
            
            data_array = hist_to_num(get_hist(data_hist, control_cat))
            mc_array = hist_to_num(get_hist(mc_hist, control_cat))
            ff_num_array = hist_to_num(get_hist(mc_hist, dr_num))
            ff_den_array = hist_to_num(get_hist(mc_hist, dr_den))
            
            qcd_array = data_array - mc_array
            qcd_values = qcd_array()
            qcd_values[qcd_array() < 0] = 0
            qcd_variances = qcd_array(sn.UP, sn.ALL, unc=True)**2
            
            cat_idx = find_by_id(qcd_hist, signal_category.id)
            qcd_hist.view().value[cat_idx, ...] = qcd_values
            qcd_hist.view().variance[cat_idx, ...] = qcd_variances     
        
        return hists
    
    
    def ff_method(task, hists):
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

        sig_reg = task.categories[0]
        
        def get_ar_data_hist(config: od.Config, hists: dict, region: str)-> hist.Hist :
            h_reg = hists[region]
            data_hists = [h for p, h in h_reg.items() if p.is_data]
            data_hist = sum(data_hists[1:], data_hists[0].copy())
            cat_id = config.get_category(region).id
            return data_hist[{"category": hist.loc(cat_id)}]
        
        def calc_wj_yields(config: od.Config, hists: dict, region: str)-> hist.Hist :
            h_reg = hists[region]
            data_hists = [h for p, h in h_reg.items() if p.is_data]
            data_hist = sum(data_hists[1:], data_hists[0].copy())
            cat_id = config.get_category(region).id
            data_ar = data_hist[{"category": hist.loc(cat_id)}]
            
            wj_hists = [h for p, h in h_reg.items() if 'wj' in p.name]
            wj_hist = sum(wj_hists[1:], wj_hists[0].copy())
            wj_ar = wj_hist[{"category": hist.loc(cat_id)}]
            
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
        for cat_idx, c in enumerate (hists_sr[jet_fakes].axes['category']):
            if c == cat_id: break
        hists_sr[jet_fakes].view().value[cat_idx,...] = fakes_val
        hists_sr[jet_fakes].view().variance[cat_idx,...] = fakes_err
        
        del hists_sr[wj_proc[0]]   
        return hists_sr
    

    config.x.hist_hooks = {
        "qcd_hook"  : qcd_estimation,
        "ff_method" : ff_method,
    }
