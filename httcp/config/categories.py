# coding: utf-8

"""
Definition of categories.
"""

import order as od

from columnflow.config_util import add_category
from columnflow.util import DotDict
from columnflow.util import maybe_import
np = maybe_import("numpy")


def add_categories(config: od.Config,
                   channel = None) -> None:
    
    
    def add_base_categories(config, channel, category_map, base_selection=[]):
        base_cat = config.get_category('_'.join(('cat',channel)))      
        for i, (cat_name, cat) in enumerate(category_map.items()):
            kwargs = {
                'name'      : '_'.join((base_cat.name, cat_name)),
                'selection' : base_selection + cat.selection,
                'id'        : 100*(i+1)+ base_cat.id,
                'label'     :' '.join((base_cat.label.split(' ')[0], cat.label if 'label' in cat else cat_name))
            }
            if 'aux' in cat.keys():
                kwargs['aux'] = {}
                for (aux_spec, aux_content) in cat.aux.items():
                    kwargs['aux'][aux_spec] = {key: '_'.join((base_cat.name, val)) for (key,val) in aux_content.items()}
                    
            add_category(config, **kwargs)
            
    
    def add_child_categories(config, parent_categories, child_category_map):
        max_child_id = int(np.max(config.categories.ids())/10000)
        for cat_name in parent_categories:
            if cat_name in ['incl', 'cat_mutau', 'cat_etau']: continue
            parent_cat = config.get_category(cat_name)
            #in case of add_child_categories function was applied previosly its needed to provide unique ids, 
            # max_child_id will be added to iterator
            for i, (child_name, child_cat) in enumerate(child_category_map.items()):
                kwargs = {
                    'name'      : '__'.join((parent_cat.name, child_name)),
                    'selection' : parent_cat.selection + child_cat.selection,
                    'id'        : 10000*(i+1+max_child_id)+ parent_cat.id,
                    'label'     : ' '.join((parent_cat.label, child_cat.label if 'label' in child_cat else child_name))
                }
                if parent_cat.aux:
                    for (aux_spec, aux_content) in parent_cat.aux:
                        if aux_spec == 'signal_reg':
                            kwargs['aux']['signal_reg'] = True
                        else:
                            kwargs['aux'][aux_spec] = aux_content

                add_category(config, **kwargs)
    
    """
    Adds all categories to a *config*.
    ids from 1 to 9 are reserved for channels
    """
    
    add_category(
        config,
        name="incl",
        id=1,
        selection=["cat_incl"],
        label="inclusive",
    )
    ##############################################
    ### Main categories for the three channels ###
    ##############################################
    
    if channel=='mutau':
        add_category(
            config,
            name="cat_mutau",
            id=2,
            selection=["cat_mutau"],
            label=r"$\mu\tau$ inclusive",
        )
        
    elif channel=='etau':
        add_category(
            config,
            name="cat_etau",
            id=3,
            selection=["cat_etau"],
            label=r"$e\tau$ inclusive",
        )
    
    #Define initial category map with selections and call the function
    #Don't change this part: it is important for fake factor method    
    base_selection = [f'cat_{channel}']
    
    category_map  = DotDict.wrap({
        "sr"            : { 'selection' : ["mt_cut", "deep_tau_wp", "lep_iso", "os_charge"],
                            'label'     : "signal region",
                            'aux'       : {
                                           #qcd estimation categories
                                           'abcd_regs' : {
                                               'ar'    : 'abcd_ar',
                                               'dr_num': 'abcd_dr_num',
                                               'dr_den':  'abcd_dr_den',
                                               },
                                           #fake factor categories
                                           'ff_regs': {
                                               "ar_wj"      : "ar_wj",
                                               "dr_num_wj"  : "dr_num_wj",
                                               "dr_den_wj"  : "dr_den_wj",
                                               "ar_qcd"      : "ar_qcd",
                                               "dr_num_qcd"  : "dr_num_qcd",
                                               "dr_den_qcd"  : "dr_den_qcd",
                                               "ar_yields"   : "ar_yields",
                                           },},},
        
        #categories for QCD estimation via classic ABCD method 
        "abcd_ar"       : { 'selection' : ["mt_cut", "deep_tau_wp", "lep_iso", "ss_charge"], 'label' : "same sign region"},
        "abcd_dr_num"   : { 'selection' : ["mt_cut", "deep_tau_wp", "lep_inv_iso", "os_charge"]},
        "abcd_dr_den"   : { 'selection' : ["mt_cut", "deep_tau_wp", "lep_inv_iso", "ss_charge"]},
        
        #categories for jet fakes estimation via classic Fake Factor method 
        "ar_wj"         : {'selection' : ["mt_cut",      "deep_tau_inv_wp",  "lep_iso", "os_charge"],},
        "dr_num_wj"     : {'selection' : ["mt_inv_cut",  "deep_tau_wp",      "lep_iso", "tau_no_fakes", "os_charge"],},
        "dr_den_wj"     : {'selection' : ["mt_inv_cut",  "deep_tau_inv_wp",  "lep_iso", "tau_no_fakes", "os_charge"],},
        
        "ar_qcd"        : {'selection' : ["mt_cut",      "deep_tau_inv_wp",  "lep_iso", "os_charge"],},
        "dr_num_qcd"    : {'selection' : ["mt_cut",      "deep_tau_wp",      "lep_inv_iso", "tau_no_fakes", "os_charge"],},
        "dr_den_qcd"    : {'selection' : ["mt_cut",      "deep_tau_inv_wp",  "lep_inv_iso", "tau_no_fakes", "os_charge"],},
        
        "ar_yields"     : {'selection' : ["mt_cut",      "deep_tau_inv_wp",  "lep_iso", "os_charge"],},
        #signal region without mt cut
        "sr_no_mt"               : { 'selection'   : ["deep_tau_wp",  "lep_iso","os_charge"],
                                 'label'        : "no mT",
                                 'aux'          : {
                                                   #qcd estimation categories
                                                        'abcd_regs' : {
                                                            'ar': 'no_mt_abcd_ar',
                                                            'dr_num'    : 'no_mt_abcd_dr_num',
                                                            'dr_den'    : 'no_mt_abcd_dr_den',
                                                        },},},
        #categories for QCD estimation via classic ABCD method 
        "no_mt_abcd_ar"       : { 'selection' : ["deep_tau_wp", "lep_iso", "ss_charge"], 'label' : "no mT, same sign region"},
        "no_mt_abcd_dr_num"   : { 'selection' : ["deep_tau_wp", "lep_inv_iso", "os_charge"]},
        "no_mt_abcd_dr_den"   : { 'selection' : ["deep_tau_wp", "lep_inv_iso", "ss_charge"]},
        })
    
    add_base_categories(config, channel, category_map, base_selection)
    #from IPython import embed; embed()
    #Add child categories to base categories
    child_category_map  = DotDict.wrap({
        "tau_barrel"    : {'selection' : ["tau_barrel"],
                            'label'     : r"$\eta_{\tau} \leq 1.2$",},
        "tau_endcap"    : {'selection' : ["tau_endcap"],
                            'label'     : r"$\eta_{\tau} > 1.2$",},
        })
    # add_child_categories(config,
    #                      parent_categories=config.categories.names(),
    #                      child_category_map=child_category_map)
