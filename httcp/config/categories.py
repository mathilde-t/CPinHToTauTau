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
            id=100 + mutau.id,
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
            id=150 + mutau.id,
            selection=["cat_mutau"  ,
                       "ss_charge"  ,
                       "mt_cut"     ,
                       "deep_tau_wp",
                       "b_veto"     ,
                       ],
            label=r"$\mu\tau$ control region",
        )
        mutau_signal_reg_no_mt = add_category(
            config,
            name="mutau_signal_reg_no_mt",
            id=200 + mutau.id,
            selection=["cat_mutau"  ,
                       "os_charge"  ,
                       "deep_tau_wp",
                       "b_veto"     ,
                       ],
            label=r"$\mu\tau$ no mt",
            aux={'control_reg': "mutau_control_reg_no_mt"} 
        )
        mutau_control_reg_no_mt = add_category(
            config,
            name="mutau_control_reg_no_mt",
            id=250 + mutau.id,
            selection=["cat_mutau"  ,
                       "ss_charge"  ,
                       "deep_tau_wp",
                       "b_veto"     ,
                       ],
            label=r"$\mu\tau$ control region no mt",
        )
        mutau_signal_reg_endcap_tau = add_category(
            config,
            name="mutau_signal_reg_endcap_tau",
            id=300 + mutau.id,
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
            id=350 + mutau.id,
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
            id=400 + mutau.id,
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
            id=450 + mutau.id,
            selection=["cat_mutau"  ,
                       "ss_charge"  ,
                       "mt_cut"     ,
                       "deep_tau_wp",
                       "b_veto"     ,
                       "tau_barrel" ,
                       ],
            label=r"$\mu\tau$ control region\n $\eta_{\tau} \leq 1.2$",
        )

        mutau_signal_reg_endcap_tau_no_mt  = add_category(
            config,
            name="mutau_signal_reg_endcap_tau_no_mt",
            id=500 + mutau.id,
            selection=["cat_mutau"  ,
                       "os_charge"  ,
                       "deep_tau_wp",
                       "b_veto"     ,
                       "tau_endcap" ,
                       ],
            label=r"$\mu\tau$ no mt \n $\eta_{\tau} > 1.2$",
            aux={'control_reg': "mutau_control_reg_endcap_tau_no_mt"}
        )
        mutau_control_reg_endcap_tau_no_mt = add_category(
            config,
            name="mutau_control_reg_endcap_tau_no_mt",
            id=550 + mutau.id,
            selection=["cat_mutau"  ,
                       "ss_charge"  ,
                       "deep_tau_wp",
                       "b_veto"     ,
                       "tau_endcap" ,
                       ],
            label=r"$\mu\tau$ no mt control region\n $\eta_{\tau} > 1.2$",
        )

        mutau_signal_reg_barrel_tau_no_mt = add_category(
            config,
            name="mutau_signal_reg_barrel_tau_no_mt",
            id=600 + mutau.id,
            selection=["cat_mutau"  ,
                       "os_charge"  ,
                       "deep_tau_wp",
                       "b_veto"     ,
                       "tau_barrel" ,
                       ],
            label=r"$\mu\tau$ no mt \n $\eta_{\tau} \leq 1.2$",
            aux={'control_reg': "mutau_control_reg_barrel_tau_no_mt"}
        )

        mutau_control_reg_barrel_tau_no_mt = add_category(
            config,
            name="mutau_control_reg_barrel_tau_no_mt",
            id=650 + mutau.id,
            selection=["cat_mutau"  ,
                       "ss_charge"  ,
                       "deep_tau_wp",
                       "b_veto"     ,
                       "tau_barrel" ,
                       ],
            label=r"$\mu\tau$ no mt control region\n $\eta_{\tau} \leq 1.2$",
        )

        mutau_signal_reg_no_mt_bveto_wp_mtt = add_category(
            config,
            name="mutau_signal_reg_no_mt_bveto_wp_mtt",
            id=700 + mutau.id,
            selection=["cat_mutau"  ,
                       "os_charge"  ,
                       "deep_tau_wp_mtt",
                      ],
            label=r"$\mu\tau$ signal region no mt, mtt wp, no b veto",
            aux={'control_reg': "mutau_contol_reg_no_mt_bveto_wp_mtt"}
        )

        mutau_contol_reg_no_mt_bveto_wp_mtt = add_category(
            config,
            name="mutau_contol_reg_no_mt_bveto_wp_mtt",
            id=750 + mutau.id,
            selection=["cat_mutau"  ,
                       "ss_charge"  ,
                       "deep_tau_wp_mtt",
                      ],
            label=r"$\mu\tau$ control region no mt, mtt wp, no b veto",
        )


        mutau_signal_reg_no_mt_bveto_wp_mtt = add_category(
            config,
            name="mutau_signal_reg_no_mt_bveto",
            id=800 + mutau.id,
            selection=["cat_mutau"  ,
                       "os_charge"  ,
                       "deep_tau_wp",
                      ],
            label=r"$\mu\tau$ signal region no mt, no b veto",
            aux={'control_reg': "mutau_contol_reg_no_mt_bveto"}
        )

        mutau_contol_reg_no_mt_bveto = add_category(
            config,
            name="mutau_contol_reg_no_mt_bveto",
            id=850 + mutau.id,
            selection=["cat_mutau"  ,
                       "ss_charge"  ,
                       "deep_tau_wp",
                      ],
            label=r"$\mu\tau$ control region no mt, no b veto",
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
            id=100 + etau.id,
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
            id=150 + etau.id,
            selection=["cat_etau"    ,
                       "ss_charge"   ,
                       "mt_cut"      ,
                       "deep_tau_wp" ,
                       "b_veto"      ,
                       ],
            label=r"$e\tau$ control region",
        )
        
        etau_signal_reg_no_mt = add_category(
            config,
            name="etau_signal_reg_no_mt",
            id=200 + etau.id,
            selection=["cat_etau"   ,
                       "os_charge"  ,
                       "deep_tau_wp",
                       "b_veto"     ,
                       ],
            label=r"$e\tau$ no mt",
            aux={'control_reg': "etau_control_reg_no_mt"}
        )
        
        etau_control_reg_no_mt = add_category(
            config,
            name="etau_control_reg_no_mt",
            id=250 + etau.id,
            selection=["cat_etau"    ,
                       "ss_charge"   ,
                       "deep_tau_wp" ,
                       "b_veto"      ,
                       ],
            label=r"$e\tau$ control region no mt",
        )
        
        etau_signal_reg_endcap_tau = add_category(
            config,
            name="etau_signal_reg_endcap_tau",
            id=300 + etau.id,
            selection=["cat_etau"   ,
                       "os_charge"  ,
                       "mt_cut"     ,
                       "deep_tau_wp",
                       "b_veto"     ,
                       "tau_endcap" ,
                       ],
            label=r"$e\tau$ signal region\n $\eta_{\tau} > 1.2$",
            aux={'control_reg': "etau_control_reg_endcap_tau"}
        )
        etau_control_reg_endcap_tau = add_category(
            config,
            name="etau_control_reg_endcap_tau",
            id=350 + etau.id,
            selection=["cat_etau"   ,
                       "ss_charge"  ,
                       "mt_cut"     ,
                       "deep_tau_wp",
                       "b_veto"     ,
                       "tau_endcap" ,
                       ],
            label=r"$e\tau$ control region\n $\eta_{\tau} > 1.2$",
        )
        
        etau_signal_reg_barrel_tau = add_category(
            config,
            name="etau_signal_reg_barrel_tau",
            id=400 + etau.id,
            selection=["cat_etau"  ,
                        "os_charge"  ,
                        "mt_cut"     ,
                        "deep_tau_wp",
                        "b_veto"     ,
                        "tau_barrel" ,
                        ],
            label=r"$e\tau$ signal region\n $\eta_{\tau} \leq 1.2$",
            aux={'control_reg': "etau_control_reg_barrel_tau"}
        )
        
        etau_control_reg_barrel_tau = add_category(
            config,
            name="etau_control_reg_barrel_tau",
            id=450 + etau.id,
            selection=["cat_etau"  ,
                        "ss_charge"  ,
                        "mt_cut"     ,
                        "deep_tau_wp",
                        "b_veto"     ,
                        "tau_barrel" ,
                        ],
            label=r"$e\tau$ control region\n $\eta_{\tau} \leq 1.2$",
        )

        etau_signal_reg_endcap_tau_no_mt = add_category(
            config,
            name="etau_signal_reg_endcap_tau_no_mt",
            id=500 + etau.id,
            selection=["cat_etau"   ,
                       "os_charge"  ,
                       "deep_tau_wp",
                       "b_veto"     ,
                       "tau_endcap" ,
                       ],
            label=r"$e\tau$ no mt signal region\n $\eta_{\tau} > 1.2$",
            aux={'control_reg': "etau_control_reg_endcap_tau_no_mt"}
        )
        etau_control_reg_endcap_tau_no_mt = add_category(
            config,
            name="etau_control_reg_endcap_tau_no_mt",
            id=550 + etau.id,
            selection=["cat_etau"   ,
                       "ss_charge"  ,
                       "deep_tau_wp",
                       "b_veto"     ,
                       "tau_endcap" ,
                       ],
            label=r"$e\tau$ no mt control region\n $\eta_{\tau} > 1.2$",
        )

        etau_signal_reg_barrel_tau_no_mt = add_category(
            config,
            name="etau_signal_reg_barrel_tau_no_mt",
            id=600 + etau.id,
            selection=["cat_etau"  ,
                        "os_charge"  ,
                        "deep_tau_wp",
                        "b_veto"     ,
                        "tau_barrel" ,
                        ],
            label=r"$e\tau$ no mt signal region\n $\eta_{\tau} \leq 1.2$",
            aux={'control_reg': "etau_control_reg_barrel_tau_no_mt"}
        )

        etau_control_reg_barrel_tau_no_mt = add_category(
            config,
            name="etau_control_reg_barrel_tau_no_mt",
            id=650 + etau.id,
            selection=["cat_etau"  ,
                        "ss_charge"  ,
                        "deep_tau_wp",
                        "b_veto"     ,
                        "tau_barrel" ,
                        ],
            label=r"$e\tau$ no mt control region\n $\eta_{\tau} \leq 1.2$",
        )


        etau_signal_reg_no_mt_bveto_wp_mtt = add_category(
            config,
            name="etau_signal_reg_no_mt_bveto_wp_mtt",
            id=700 + etau.id,
            selection=["cat_etau"  ,
                       "os_charge"  ,
                       "deep_tau_wp_mtt",
                      ],
            label=r"$e\tau$ signal region no mt, mtt wp, no b veto",
            aux={'control_reg': "etau_contol_reg_no_mt_bveto_wp_mtt"}
        )

        etau_contol_reg_no_mt_bveto_wp_mtt = add_category(
            config,
            name="etau_contol_reg_no_mt_bveto_wp_mtt",
            id=750 + etau.id,
            selection=["cat_etau"  ,
                       "ss_charge"  ,
                       "deep_tau_wp_mtt",
                      ],
            label=r"$e\tau$ control region no mt, mtt wp, no b veto",
        )


        etau_signal_reg_no_mt_bveto_wp_mtt = add_category(
            config,
            name="etau_signal_reg_no_mt_bveto",
            id=800 + etau.id,
            selection=["cat_etau"  ,
                       "os_charge"  ,
                       "deep_tau_wp",
                      ],
            label=r"$e\tau$ signal region no mt, no b veto",
            aux={'control_reg': "etau_contol_reg_no_mt_bveto"}
        )

        etau_contol_reg_no_mt_bveto = add_category(
            config,
            name="etau_contol_reg_no_mt_bveto",
            id=850 + etau.id,
            selection=["cat_etau"  ,
                       "ss_charge"  ,
                       "deep_tau_wp",
                      ],
            label=r"$e\tau$ control region no mt, no b veto",
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


    

  
