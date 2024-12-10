# coding: utf-8

"""
Definition of categories.
"""

import order as od

from columnflow.config_util import add_category


def add_categories(config: od.Config,
                   channel = None) -> None:
    """
    Adds all categories to a *config*.
    ids from 1 to 9 are reserved for channels
    ids from 10 to 990 are reserved for helper categories like mt cut or same/opposite sign selection
    ids from 1000 to 990000 are reserved for actual signal, application and determination regions
    """
    
    
    add_category(
        config,
        name="incl",
        id=1,
        selection="cat_incl",
        label="inclusive",
    )
     
    ##############################################
    ### Main categories for the three channels ###
    ##############################################
    if channel=='mutau':
        mutau = add_category(
            config,
            name="cat_mutau",
            id=2,
            selection="cat_mutau",
            label=r"$\mu\tau$ inclusive",
        )
        #################################
        ### mu-tau channel categories ###
        #################################
        
        mutau_signal_reg = add_category(
            config,
            name="mutau_signal_reg",
            id=2000 + mutau.id,
            selection=["cat_mutau"  ,
                       "os_charge"  ,
                       "mt_cut"     ,
                       "deep_tau_wp",
                       "b_veto"     ,
                      ],
            label=r"$\mu\tau$ signal region",
            aux={'control_reg': "mutau_control_reg"}
        )
        
        mutau_control_reg = add_category(
            config,
            name="mutau_control_reg",
            id=5000 + mutau.id,
            selection=["cat_mutau"  ,
                       "ss_charge"  ,
                       "mt_cut"     ,
                       "deep_tau_wp",
                       "b_veto"     ,
                       ],
            label=r"$\mu\tau$ control region",
        )
        mutau_no_mt = add_category(
            config,
            name="mutau_no_mt",
            id=3000 + mutau.id,
            selection=["cat_mutau"  ,
                       "os_charge"  ,
                       "deep_tau_wp",
                       "b_veto"     ,
                       ],
            label=r"$mu\tau$ no mt",
        )
        mutau_signal_reg_endcap_tau = add_category(
            config,
            name="mutau_signal_reg_endcap_tau",
            id=6000 + mutau.id,
            selection=["cat_mutau"  ,
                       "os_charge"  ,
                       "mt_cut"     ,
                       "deep_tau_wp",
                       "b_veto"     ,
                       "tau_endcap" ,
                       ],
            label=r"$\mu\tau$ signal region\n $\eta_{\tau} > 1.2$",
            aux={'control_reg': "mutau_control_reg_endcap_tau"}
        )
        mutau_control_reg_endcap_tau = add_category(
            config,
            name="mutau_control_reg_endcap_tau",
            id=7000 + mutau.id,
            selection=["cat_mutau"  ,
                       "ss_charge"  ,
                       "mt_cut"     ,
                       "deep_tau_wp",
                       "b_veto"     ,
                       "tau_endcap" ,
                       ],
            label=r"$\mu\tau$ control region\n $\eta_{\tau} > 1.2$",
        )
        
        mutau_signal_reg_barrel_tau = add_category(
            config,
            name="mutau_signal_reg_barrel_tau",
            id=8000 + mutau.id,
            selection=["cat_mutau"  ,
                       "os_charge"  ,
                       "mt_cut"     ,
                       "deep_tau_wp",
                       "b_veto"     ,
                       "tau_barrel" ,
                       ],
            label=r"$\mu\tau$ signal region\n $\eta_{\tau} \leq 1.2$",
            aux={'control_reg': "mutau_control_reg_barrel_tau"}
        )
        
        mutau_control_reg_barrel_tau = add_category(
            config,
            name="mutau_control_reg_barrel_tau",
            id=9000 + mutau.id,
            selection=["cat_mutau"  ,
                       "ss_charge"  ,
                       "mt_cut"     ,
                       "deep_tau_wp",
                       "b_veto"     ,
                       "tau_barrel" ,
                       ],
            label=r"$\mu\tau$ control region\n $\eta_{\tau} \leq 1.2$",
        )

    elif channel=='etau':
        etau = add_category(
            config,
            name="cat_etau",
            id=3,
            selection="cat_etau",
            label=r"$e\tau$ inclusive",
        )
        ################################
        ### e-tau channel categories ###
        ################################
        etau_signal_reg = add_category(
            config,
            name="etau_signal_reg",
            id=3000 + etau.id,
            selection=["cat_etau"  ,
                       "os_charge"  ,
                       "mt_cut"     ,
                       "deep_tau_wp",
                       "b_veto"     ,
                       ],
            label=r"$e\tau$ signal region",
            aux={'control_reg': "etau_control_reg"}
        )
        etau_control_reg = add_category(
            config,
            name="etau_control_reg",
            id=6000 + etau.id,
            selection=["cat_etau"    ,
                       "ss_charge"   ,
                       "mt_cut"      ,
                       "deep_tau_wp" ,
                       "b_veto"      ,
                       ],
            label=r"$e\tau$ control region",
        )
        
        etau_no_mt = add_category(
            config,
            name="etau_no_mt",
            id=2000 + etau.id,
            selection=["cat_etau"  ,
                       "os_charge"  ,
                       "deep_tau_wp",
                       "b_veto"     ,
                       ],
            label=r"$e\tau$ no mt",
        )
    
    elif channel=='tautau':
        tautau = add_category(
            config,
            name="cat_tautau",
            id=4,
            selection="cat_tautau",
            label=r"$\tau\tau$ inclusive",
        )

    else:
        raise ValueError(f"attempt to process more than one channel: {ch_str}")


    

  
