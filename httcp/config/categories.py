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
    mutau = add_category(
        config,
        name="cat_mutau",
        id=2,
        selection="cat_mutau",
        label=r"$\mu\tau$ inclusive",
    )
    
    # etau = add_category(
    #     config,
    #     name="cat_etau",
    #     id=3,
    #     selection="cat_etau",
    #     label=r"$e\tau$ inclusive",
    # )
    
    # tautau = add_category(
    #     config,
    #     name="cat_tautau",
    #     id=4,
    #     selection="cat_tautau",
    #     label=r"$\tau\tau$ inclusive",
    # )
    ##############################################
    ### Main categories for the three channels ###
    ##############################################
    
    #############################################
    ### Categories for the Fake Factor method ###
    #############################################
    
    # os_charge = add_category(
    #     config,
    #     name="os_charge",
    #     id=10,
    #     selection="os_charge",
    #     label=r"$q_1\cdot q_2 < 0$",
    # )
    
    # ss_charge = add_category(
    #     config,
    #     name="ss_charge",
    #     id=20,
    #     selection="ss_charge",
    #     label=r"$q_1\cdot q_2 > 0$",
    # )
    
    # mT_cut = add_category(
    #     config,
    #     name="mt_cut",
    #     id=30,
    #     selection="mt_cut",
    #     label=r"$mT < 50$ GeV",
    # )
    
    # mT_inv_cut = add_category(
    #     config,
    #     name="mt_inv_cut",
    #     id=40,
    #     selection="mt_inv_cut",
    #     label=r"$mT \geq 50$ GeV",
    # )
    
    # deep_tau_wp = add_category(
    #     config,
    #     name="deep_tau_wp",
    #     id=50,
    #     selection="deep_tau_wp",
    #     label="DeepTau wp",
    # )
    # deep_tau_inv_wp = add_category(
    #     config,
    #     name="deep_tau_inv_wp",
    #     id=60,
    #     selection="deep_tau_inv_wp",
    #     label="inv DeepTau",
    # )
    # b_veto = add_category(
    #     config,
    #     name="b_veto",
    #     id=70,
    #     selection="b_veto",
    #     label="b veto",
    # )
    # b_veto_inv = add_category(
    #     config,
    #     name="b_veto_inv",
    #     id=80,
    #     selection="b_veto_inv",
    #     label="inv b veto",
    # )
    #############################################
    ### Categories for the Fake Factor method ###
    #############################################
    
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

    elif channel=='etau':
        etau = add_category(
            config,
            name="cat_etau",
            id=3,
            selection="cat_etau",
            label=r"$e\tau$ inclusive",
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
    ##############################################
    ### Main categories for the three channels ###
    ##############################################
    
    #################################
    ### mu-tau channel categories ###
    #################################
    
    mutau_signal_reg = add_category(
        config,
        name="mutau_signal_reg",
        id=2000 + mutau.id,
        selection=["cat_mutau",
                   "os_charge",
                   "mt_cut",
                   "deep_tau_wp",
                   "b_veto",
                   ],
        label=r"$\mu\tau$ signal region",
        aux={'control_reg': "mutau_control_reg"}
    )
    
    mutau_control_reg = add_category(
        config,
        name="mutau_control_reg",
        id=5000 + mutau.id,
        selection=["cat_mutau",
                   "ss_charge",
                   "mt_cut",
                   "deep_tau_wp",
                   "b_veto",
                   ],
        label=r"$\mu\tau$ control region",
    )
    
    mutau_signal_reg_endcap_tau = add_category(
        config,
        name="mutau_signal_reg_endcap_tau",
        id=6000 + mutau.id,
        selection=["cat_mutau",
                   "os_charge",
                   "mt_cut",
                   "deep_tau_wp",
                   "b_veto",
                   "tau_endcap",
                   ],
        label=r"$\mu\tau$ signal region\n $\eta_{\tau} > 1.2$",
        aux={'control_reg': "mutau_control_reg_endcap_tau"}
    )
    mutau_control_reg_endcap_tau = add_category(
        config,
        name="mutau_control_reg_endcap_tau",
        id=7000 + mutau.id,
        selection=["cat_mutau",
                   "ss_charge",
                   "mt_cut",
                   "deep_tau_wp",
                   "b_veto",
                   "tau_endcap",
                   ],
        label=r"$\mu\tau$ control region\n $\eta_{\tau} > 1.2$",
    )
    
    mutau_signal_reg_barrel_tau = add_category(
        config,
        name="mutau_signal_reg_barrel_tau",
        id=8000 + mutau.id,
        selection=["cat_mutau",
                   "os_charge",
                   "mt_cut",
                   "deep_tau_wp",
                   "b_veto",
                   "tau_barrel",
                   ],
        label=r"$\mu\tau$ signal region\n $\eta_{\tau} \leq 1.2$",
        aux={'control_reg': "mutau_control_reg_barrel_tau"}
    )
    
    mutau_control_reg_barrel_tau = add_category(
        config,
        name="mutau_control_reg_barrel_tau",
        id=9000 + mutau.id,
        selection=["cat_mutau",
                   "ss_charge",
                   "mt_cut",
                   "deep_tau_wp",
                   "b_veto",
                   "tau_barrel",
                   ],
        label=r"$\mu\tau$ control region\n $\eta_{\tau} \leq 1.2$",
    )
    
    #################################
    ### mu-tau channel categories ###
    #################################
    
    ################################
    ### e-tau channel categories ###
    ################################

    elif channel=='etau':
        etau_signal_reg = add_category(
            config,
            name="etau_signal_reg",
            id=3000 + etau.id,
            selection=[etau.selection,
                       os_charge.selection,
                       mT_cut.selection,
                       deep_tau_wp.selection,
                       b_veto.selection],
            label=r"$e\tau$ signal region",
        )
        etau_ff_control_reg = add_category(
            config,
            name="etau_ff_control_reg",
            id=5000 + etau.id,
            selection=[etau.selection,
                       ss_charge.selection,
                       mT_cut.selection,
                       deep_tau_wp.selection,
                       b_veto.selection],
            label=r"$e\tau$ control region",
        )
        etau_inv_deeptau = add_category(
            config,
            name="etau_inv_deeptau",
            id=4000 + etau.id,
            selection=[etau.selection,
                       os_charge.selection,
                       mT_cut.selection,
                       b_veto.selection,
                       deep_tau_inv_wp.selection],
            label=r"$e\tau$ inv DeepTau",
        )
    
    # etau_signal_reg = add_category(
    #     config,
    #     name="etau_signal_reg",
    #     id=3000 + etau.id,
    #     selection=[etau.selection,
    #                os_charge.selection,
    #                mT_cut.selection,
    #                deep_tau_wp.selection,
    #                b_veto.selection],
    #     label=r"$e\tau$ signal region",
    # )
    # etau_ff_control_reg = add_category(
    #     config,
    #     name="etau_ff_control_reg",
    #     id=5000 + etau.id,
    #     selection=[etau.selection,
    #                ss_charge.selection,
    #                mT_cut.selection,
    #                deep_tau_wp.selection,
    #                b_veto.selection],
    #     label=r"$e\tau$ control region",
    # )
    # etau_inv_deeptau = add_category(
    #     config,
    #     name="etau_inv_deeptau",
    #     id=4000 + etau.id,
    #     selection=[etau.selection,
    #                os_charge.selection,
    #                mT_cut.selection,
    #                b_veto.selection,
    #                deep_tau_inv_wp.selection],
    #     label=r"$e\tau$ inv DeepTau",
    # )
    
    # etau_no_mt = add_category(
    #     config,
    #     name="etau_no_mt",
    #     id=2000 + etau.id,
    #     selection=[etau.selection,
    #                os_charge.selection,
    #                b_veto.selection,
    #                deep_tau_wp.selection],
    #     label=r"$e\tau$ no mt",
    # )
    
    ################################
    ### e-tau channel categories ###
    ################################
