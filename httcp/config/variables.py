### This config is used for listing the variables used in the analysis ###

from columnflow.config_util import add_category

import order as od

from columnflow.columnar_util import EMPTY_FLOAT
from columnflow.util import DotDict
from columnflow.columnar_util import ColumnCollection

from columnflow.util import maybe_import
np = maybe_import("numpy")

def keep_columns(cfg: od.Config) -> None:
    # columns to keep after certain steps
    cfg.x.keep_columns = DotDict.wrap({
        "cf.ReduceEvents": {
            # TauProds
            "TauProd.*",
            "GenPart.*",
            "GenZ.*",
            # general event info
            "run", "luminosityBlock", "event",
            "PV.npvs","Pileup.nTrueInt","Pileup.nPU","genWeight", "LHEWeight.originalXWGTUP", "HTXS_njets*", "LHE_Njets","weight","zpt_weight","muon_weight_nom","mc_weight","tau_weight_nom",
        } | {
            f"PuppiMET.{var}" for var in [
                "pt", "phi", "significance",
                "covXX", "covXY", "covYY", "pt_no_jec", "eta_no_jec"
            ]
        } | {
            f"MET.{var}" for var in [
                "pt", "phi", "significance",
                "covXX", "covXY", "covYY",
            ]     
        } | {
            f"Jet.{var}" for var in [
                "pt", "eta", "phi", "mass", 
                "btagDeepFlavB", "hadronFlavour", "pt_no_jec", "phi_no_jec","eta_no_jec", "mass_no_jec", "jec_no_jec_diff","number_of_jets",
            ] 
        } | {
            f"Tau.{var}" for var in [
                "pt","eta","phi","mass","dxy","dz", "charge", 
                "rawDeepTau2018v2p5VSjet","idDeepTau2018v2p5VSjet", "idDeepTau2018v2p5VSe", "idDeepTau2018v2p5VSmu", 
                "decayMode", "decayModePNet", "genPartFlav", "rawIdx",
                "pt_no_tes", "mass_no_tes", "IPx", "IPy", "IPz","ip_sig", "jetIdx"
            ] 
        } | {
            f"Muon.{var}" for var in [
                "pt","eta","phi","mass","dxy","dz", "charge",
                "decayMode", "pfRelIso04_all","mT", "rawIdx","IPx", "IPy", "IPz","ip_sig"
            ] 
        } | {
            f"Electron.{var}" for var in [
                "pt","eta","phi","mass","dxy","dz", "charge", 
                "decayMode", "pfRelIso03_all", "mT", "rawIdx", "IPx", "IPy", "IPz","ip_sig", "pt_no_scaling_smearing",
            ] 
        } | {
            f"{var}_triggerd" for var in [ #Trigger variables to have a track of a particular trigger fired
                "single_electron", "cross_electron",
                "single_muon", "cross_muon",
                "cross_tau",
            ]
        } | {
            f"matched_triggerID_{var}" for var in [
                "e", "mu", "tau",
            ]
        } | {
            f"TrigObj.{var}" for var in [
                "id", "pt", "eta", "phi", "filterBits",
            ]
        } | {
            f"TauSpinner.weight_cp_{var}" for var in [
                "0", "0_alt", "0p25", "0p25_alt", "0p375",
                "0p375_alt", "0p5", "0p5_alt", "minus0p25", "minus0p25_alt"
            ]
        } | {
            f"hcand.{var}" for var in [
                "pt","eta","phi","mass", "charge", 
                "decayMode", "rawIdx", "ip_sig", "IPx", "IPy","IPz"
            ]
        } | {
            "GenTau.*", "GenTauProd.*", "nJet",
        } | {
            f"hcandprod.{var}" for var in [
                "pt", "eta", "phi", "mass", "charge",
                "pdgId", "tauIdx"
            ]
        } | {
		"hcand_*","tau_decay_prods*"
	} | {"is_b_vetoed","channel_id"} | {ColumnCollection.ALL_FROM_SELECTOR},
        "cf.MergeSelectionMasks": {
            "normalization_weight", 
            "cutflow.*", "process_id", "category_ids",
        },
        "cf.UniteColumns": {
            "*",
        },
    })



def add_common_features(cfg: od.config) -> None:
    """
    Adds common features
    """
    cfg.add_variable(
        name="event",
        expression="event",
        binning=(1, 0.0, 1.0e9),
        x_title="Event number",
        discrete_x=True,
    )
    cfg.add_variable(
        name="N_events",
        expression="N_events",
        binning=(1, 0.0, 1.0e9),
        x_title="Event number",
        discrete_x=True,
    )
    cfg.add_variable(
        name="run",
        expression="run",
        binning=(1, 100000.0, 500000.0),
        x_title="Run number",
        discrete_x=True,
    )
    cfg.add_variable(
        name="lumi",
        expression="luminosityBlock",
        binning=(1, 0.0, 5000.0),
        x_title="Luminosity block",
        discrete_x=True,
    )


def add_lepton_features(cfg: od.Config) -> None:
    """
    Adds lepton features only , ex electron_1_pt
    """
    cfg.add_variable(
        name=f"electron_1_pt_no_scaling_smearing",
        expression="Electron.pt_no_scaling_smearing[:,0]",
        null_value=EMPTY_FLOAT,
        binning=(40, 0., 200.),
        unit="GeV",
        x_title= r" Electron $p_{T}$ no scaling or smearing",
    )
    
    for obj in ["Electron", "Muon", "Tau"]:
        for i in range(2):
            cfg.add_variable(
                name=f"{obj.lower()}_{i+1}_pt",
                expression=f"{obj}.pt[:,{i}]",
                null_value=EMPTY_FLOAT,
                binning=(40, 0., 200.),
                unit="GeV",
                x_title=obj + r" $p_{T}$",
            )
            cfg.add_variable(
                name=f"{obj.lower()}_{i+1}_phi",
                expression=f"{obj}.phi[:,{i}]",
                null_value=EMPTY_FLOAT,
                binning=(32, -3.2, 3.2),
                x_title=obj + r" $\phi$",
            )
            cfg.add_variable(
                name=f"{obj.lower()}_{i+1}_eta",
                expression=f"{obj}.eta[:,{i}]",
                null_value=EMPTY_FLOAT,
                binning=(25, -2.5, 2.5),
                x_title=obj + r" $\eta$",
            )
        cfg.add_variable(
            name=f"{obj.lower()}_ip_sig",
            expression=f"{obj}.ip_sig",
            null_value=EMPTY_FLOAT,
            binning=(40, 0.0, 10),
            unit="",
            x_title=obj + r"$\frac{|IP|}{\sigma(IP)}$",
        )


def add_jet_features(cfg: od.Config) -> None:
    """
    Adds jet features only
    """
    cfg.add_variable(
        name="n_jet",
        expression="nJet",
        binning=(11, -0.5, 10.5),
        x_title="Number of jets",
        discrete_x=True,
    )
    cfg.add_variable(
        name="jets_pt_no_jec",
        expression="Jet.pt_no_jec",
        binning=(40, 0.0, 400.0),
        unit="GeV",
        x_title=r"$Uncorrected p_{T} of all jets$",
    )      
    cfg.add_variable(
        name="jets_pt",
        expression="Jet.pt",
        binning=(40, 0.0, 400.0),
        unit="GeV",
        x_title=r"$p_{T} of all jets$",
    )
    cfg.add_variable(
        name="Number_of_Jets",
        expression="Jet.number_of_jets",
        binning=(11, -0.5, 10.5),
        x_title=r"Number of Jets",
        discrete_x=True,
    )
    for i in range(2):
        cfg.add_variable(
            name=f"jet_{i+1}_pt",
            expression=f"Jet.pt[:,{i}]",
            null_value=EMPTY_FLOAT,
            binning=(40, 0.0, 400.0),
            unit="GeV",
            x_title=r"Jet $p_{T}$",
        )
        cfg.add_variable(
            name=f"jet_{i+1}_eta",
            expression=f"Jet.eta[:,{i}]",
            null_value=EMPTY_FLOAT,
            binning=(30, -3., 3.),
            x_title=r"Jet $\eta$",
        )
        cfg.add_variable(
            name=f"jet_{i+1}_phi",
            expression=f"Jet.phi[:,{i}]",
            null_value=EMPTY_FLOAT,
            binning=(32, -3.2, 3.2),
            x_title=r"Jet $\varphi$",
        )
        cfg.add_variable(
            name=f"jet_{i+1}_pt_no_jec",
            expression=f"Jet.pt_no_jec[:,{i}]",
            null_value=EMPTY_FLOAT,
            binning=(40, 0.0, 400.0),
            unit="GeV",
            x_title=r"Jet $p_{T} no jec $",
        )
        cfg.add_variable(
            name=f"jet_{i+1}_eta_no_jec",
            expression=f"Jet.eta_no_jec[:,{i}]",
            null_value=EMPTY_FLOAT,
            binning=(30, -3.0, 3.0),
            x_title=r"Jet $\eta$ no jec ",
        )
        cfg.add_variable(
            name=f"jet_{i+1}_phi_no_jec",
            expression=f"Jet.phi_no_jec[:,{i}]",
            null_value=EMPTY_FLOAT,
            binning=(32, -3.2, 3.2),
            x_title=r"Jet $\phi$ no jec",
        )
    cfg.add_variable(
        name="jet_jec_no_jec_diff",
        expression="Jet.jec_no_jec_diff",
        null_value=EMPTY_FLOAT,
        binning=(20,-10,10),
        x_title=r"$Jet_{jec} $p_{T} - $Jet_{no jec} $p_{T}",
    )
    cfg.add_variable(
        name="ht",
        # expression=lambda events: ak.sum(events.Jet.pt, axis=1),
        expression="ht",
        binning=(40, 0.0, 800.0),
        unit="GeV",
        x_title="HT",
    )
    cfg.add_variable(
        name="jet_raw_DeepJetFlavB",
        expression="Jet.btagDeepFlavB",
        null_value=EMPTY_FLOAT,
        binning=(30, 0,1),
        x_title=r"raw DeepJetFlawB",
    )


def add_highlevel_features(cfg: od.Config) -> None:    
    """
    Adds MET and other high-level features
    """
    cfg.add_variable(
        name="met",
        expression="MET.pt",
        null_value=EMPTY_FLOAT,
        binning=(20, 0.0, 100.0),
        x_title=r"MET",
    )
    cfg.add_variable(
        name="puppi_met_pt_no_jec",
        expression="PuppiMET.pt_no_jec",
        null_value=EMPTY_FLOAT,
        binning=(50, 0,100),
        unit="GeV",
        x_title=r"Uncorrected PuppiMET $p_T$",
    )
    cfg.add_variable(
        name="puppi_met_phi_no_jec",
        expression="PuppiMET.phi_no_jec",
        null_value=EMPTY_FLOAT,
        binning=(32, -3.2,3.2),
        x_title=r"Uncorrected PuppiMET $\phi$",
    )
    cfg.add_variable(
        name="puppi_met_pt",
        expression="PuppiMET.pt",
        null_value=EMPTY_FLOAT,
        binning=(50, 0,100),
        unit="GeV",
        x_title=r"PUPPI MET $p_T$",
    )
    cfg.add_variable(
        name="puppi_met_phi",
        expression="PuppiMET.phi",
        null_value=EMPTY_FLOAT,
        binning=(32, -3.2,3.2),
        x_title=r"PUPPI MET $\phi$",
    )  
    
    
    

def add_weight_features(cfg: od.Config) -> None:
    """
    Adds weights
    """
    cfg.add_variable(
        name="mc_weight",
        expression="mc_weight",
        binning=(20, -2, 2),
        x_title="MC weight",
    )
    cfg.add_variable(
        name="pu_weight",
        expression="pu_weight",
        null_value=EMPTY_FLOAT,
        binning=(30, 0,3),
        unit="",
        x_title=r"Pileup weight",
    )
    
    cfg.add_variable(
        name="muon_weight",
        expression="muon_weight_nom",
        null_value=EMPTY_FLOAT,
        binning=(50, 0.5,1.5),
        unit="",
        x_title=r"muon weight",
    )
    
    cfg.add_variable(
        name="tau_weight",
        expression="tau_weight_nom",
        null_value=EMPTY_FLOAT,
        binning=(50, 0.5,1.5),
        unit="",
        x_title=r"tau weight",
    )
    
    for var in ["0", "0p25", "0p375", "0p5", "minus0p25"]:
        
        angle = float(var.replace("minus","-").replace("p", "."))*180
        cfg.add_variable(
        name=f"TauSpinner_weight_cp_{var}",
        expression=f"TauSpinner.weight_cp_{var}",
        null_value=EMPTY_FLOAT,
        binning=(60, -3,3),
        unit="",
        x_title=fr"Tau spinner weight $\Delta \phi$=${angle}^{{\circ}}$",
    )
        

    

def add_cutflow_features(cfg: od.Config) -> None:
    """
    Adds cf features
    """
    cfg.add_variable(
        name="cf_jet1_pt",
        expression="cutflow.jet1_pt",
        binning=(40, 0.0, 400.0),
        unit="GeV",
        x_title=r"Jet 1 $p_{T}$",
    )


def phi_cp_variables(cfg: od.Config) -> None:
    n_bins_phi_cp = 11
    for the_ch in ['mu_pi', 'mu_rho', 'mu_a1_1pr', "rho_rho","pi_pi"]:
        spitted_str = the_ch.split('_')
        if 'a1' in the_ch: 
            title_str = "\\" + spitted_str[0] + fr" a_1, {spitted_str[2]}"
        else:
             title_str = "\\" + spitted_str[0] + "\\" + spitted_str[1]
        cfg.add_variable(
            name=f"phi_cp_{the_ch}",
            expression=f"phi_cp_{the_ch}",
            null_value=EMPTY_FLOAT,
            binning=(n_bins_phi_cp, 0, 2*np.pi),
            x_title=rf"$\varphi_{{CP}} [{title_str}]$ (rad)",
        )
        cfg.add_variable(
            name=f"phi_cp_{the_ch}_reg1",
            expression=f"phi_cp_{the_ch}_reg1",
            null_value=EMPTY_FLOAT,
            binning=(n_bins_phi_cp, 0, 2*np.pi),
            x_title=rf"$\varphi_{{CP}} [{title_str}], \alpha < \pi/4$ (rad)",
        )
        cfg.add_variable(
            name=f"phi_cp_{the_ch}_reg2",
            expression=f"phi_cp_{the_ch}_reg2",
            null_value=EMPTY_FLOAT,
            binning=(n_bins_phi_cp, 0, 2*np.pi),
            x_title=rf"$\varphi_{{CP}} [{title_str}], \alpha \geq \pi/4$ (rad)",
        )
        # 2-bin histograms
#        cfg.add_variable(
#            name=f"phi_cp_{the_ch}_2bin",
#            expression=f"phi_cp_{the_ch}_2bin",
#            null_value=EMPTY_FLOAT,
#            binning=(2, 0, 2*np.pi), 
#            x_title=rf"$\varphi_{{CP}} [{title_str}]$ (rad)",
#        )
#        cfg.add_variable(
#            name=f"phi_cp_{the_ch}_reg1_2bin",
#            expression=f"phi_cp_{the_ch}_reg1_2bin",
#            null_value=EMPTY_FLOAT,
#            binning=(2, 0, 2*np.pi),
#            x_title=rf"$\varphi_{{CP}} [{title_str}], \alpha < \pi/4$ (rad)",
#        )
#        cfg.add_variable(
#            name=f"phi_cp_{the_ch}_reg2_2bin",
#            expression=f"phi_cp_{the_ch}_reg2_2bin",
#            null_value=EMPTY_FLOAT,
#            binning=(2, 0, 2*np.pi),
#            x_title=rf"$\varphi_{{CP}} [{title_str}], \alpha \geq \pi/4$ (rad)",
#        )
#        cfg.add_variable(
#            name=f"alpha_{the_ch}",
#            expression=f"alpha_{the_ch}",
#            null_value=EMPTY_FLOAT,
#            binning=(6, 0, np.pi/2),
#            x_title=rf"$ \alpha [{title_str}] $(rad)",
#        )

def add_dilepton_features(cfg: od.Config) -> None:
    channels = cfg.channels.names()
    ch_objects = DotDict.wrap({
        'etau'   : {'lep0' : 'Electron',
                    'lep1' : 'Tau'},
        'mutau'  : {'lep0' : 'Muon',
                    'lep1' : 'Tau'},
        'tautau' : {'lep0' : 'Tau',
                    'lep1' : 'Tau'},
    })
    bin_split_factor = 4 #Used to define histograms for kinematic variables with finer binning
    for ch_str in channels:
        cfg.add_variable(
                name=f"{ch_str}_mvis",
                expression=f"hcand_{ch_str}.mass",
                null_value=EMPTY_FLOAT,
                binning=(40, 0.0, 200.0),
                unit="GeV",
                x_title=r"$m_{vis}$",
            )
        if ch_str != 'tautau':
            cfg.add_variable(
                name=f"{ch_str}_mt",
                expression=f"hcand_{ch_str}.mt",
                null_value=EMPTY_FLOAT,
                binning=(40, 0.0, 200.0),
                unit="GeV",
                x_title=r"$m_{T}$",
            )
        cfg.add_variable(
            name=f"{ch_str}_delta_r",
            expression=f"hcand_{ch_str}.delta_r",
            null_value=EMPTY_FLOAT,
            binning=(40, 0, 4),
            x_title=r"$\Delta R(\ell,\ell)$",
        )
        cfg.add_variable(
                name=f"{ch_str}_pt",
                expression=f"hcand_{ch_str}.pt",
                null_value=EMPTY_FLOAT,
                binning=(40, 0.0, 200.0),
                unit="GeV/c",
                x_title=r"$p_{T}(\ell\ell)$",
        )
        
        for lep in ['lep0','lep1']:
            if ch_str != 'tautau': lep_str = ch_objects[ch_str][lep].lower()
            else: lep_str = f'tau {lep[3:]}'
            cfg.add_variable(
                name=f"{ch_str}_{lep}_pt",
                expression=f"hcand_{ch_str}.{lep}.pt",
                null_value=EMPTY_FLOAT,
                binning=(30, 20, 80.),
                unit="GeV",
                x_title= rf"{lep_str} $p_{{T}}$",
            )
            cfg.add_variable(
                name=f"{ch_str}_{lep}_eta",
                expression=f"hcand_{ch_str}.{lep}.eta",
                null_value=EMPTY_FLOAT,
                binning=(32, -3.2, 3.2),
                x_title=rf"{lep_str} $\eta$",
            )
            cfg.add_variable(
                name=f"{ch_str}_{lep}_phi",
                expression=f"hcand_{ch_str}.{lep}.phi",
                null_value=EMPTY_FLOAT,
                binning=(32, -3.2, 3.2),
                x_title=rf"{lep_str} $\phi$",
            )
            cfg.add_variable(
                name=f"{ch_str}_{lep}_mass",
                expression=f"hcand_{ch_str}.{lep}.mass",
                null_value=EMPTY_FLOAT,
                binning=(30, 0, 3),
                unit="GeV",
                x_title=f"{lep_str} mass",
            )
            cfg.add_variable(
                name=f"{ch_str}_{lep}_decayModePNet",
                expression=f"hcand_{ch_str}.{lep}.decayModePNet",
                null_value=EMPTY_FLOAT,
                binning=(12,0,12),
                unit="",
                x_title=rf"{lep_str} PNet decay mode",
            )
            cfg.add_variable(
                name=f"{ch_str}_{lep}_decayMode",
                expression=f"hcand_{ch_str}.{lep}.decayMode",
                null_value=EMPTY_FLOAT,
                binning=(12,0,12),
                unit="",
                x_title=rf"{lep_str} HPS decay mode",
            )
            cfg.add_variable(
                name=f"{ch_str}_{lep}_ip_sig",
                expression=f"hcand_{ch_str}.{lep}.ip_sig",
                null_value=EMPTY_FLOAT,
                binning=(40, 0.0, 10),
                unit="",
                x_title= rf"{lep_str} $\frac{{|IP|}}{{\sigma(IP)}}$",
            )
            for proj in ['x','y','z']:
                cfg.add_variable(
                    name=f"{ch_str}{lep}_ip_{proj}",
                    expression=f"hcand_{ch_str}.{lep}.IP{proj}",
                    null_value=EMPTY_FLOAT,
                    binning=(30, -0.002, 0.002),
                    unit="",
                    x_title= rf"{lep_str} $IP_{proj}$",
                )
            #Variables with finer binning
            cfg.add_variable(
                name=f"{ch_str}_{lep}_pt_fine_binning",
                expression=f"hcand_{ch_str}.{lep}.pt",
                null_value=EMPTY_FLOAT,
                binning=(30*bin_split_factor, 20, 80.),
                unit="GeV",
                x_title= rf"{lep_str} $p_{{T}}$",
            )
            cfg.add_variable(
                name=f"{ch_str}_{lep}_eta_fine_binning",
                expression=f"hcand_{ch_str}.{lep}.eta",
                null_value=EMPTY_FLOAT,
                binning=(32*bin_split_factor, -3.2, 3.2),
                x_title=rf"{lep_str} $\eta$",
            )
            cfg.add_variable(
                name=f"{ch_str}_{lep}_phi_fine_binning",
                expression=f"hcand_{ch_str}.{lep}.phi",
                null_value=EMPTY_FLOAT,
                binning=(32*bin_split_factor, -3.2, 3.2),
                x_title=rf"{lep_str} $\phi$",
            )


def add_variables(cfg: od.Config) -> None:
    """
    Adds all variables to a *config*.
    """
    add_common_features(cfg)
    add_lepton_features(cfg)
    add_jet_features(cfg)
    add_highlevel_features(cfg)
    #phi_cp_variables(cfg)
    add_weight_features(cfg)
    add_cutflow_features(cfg)
    add_dilepton_features(cfg)
