# coding: utf-8

"""
Configuration of the higgs_cp analysis.
"""

import functools
import yaml
import law
import order as od
from scinum import Number

from columnflow.util import DotDict, maybe_import, dev_sandbox
from columnflow.config_util import (
    get_root_processes_from_campaign, 
    add_category, add_shift_aliases,
    verify_config_processes,get_shifts_from_sources
)

ak = maybe_import("awkward")


def add_run3_2022_preEE (ana: od.Analysis,
                      campaign: od.Campaign,
                      config_name           = None,
                      config_id             = None,
                      limit_dataset_files   = None,) -> od.Config :

    # get all root processes
    procs = get_root_processes_from_campaign(campaign)
    
    # create a config by passing the campaign, so id and name will be identical
    cfg = ana.add_config(campaign,
                        name  = config_name,
                        id    = config_id)

    # gather campaign data
    cfg.x.year = campaign.x.year
    
    # validations
    if campaign.x.year == 2022:
        assert campaign.x.tag in ["preEE", "postEE"]
    # gather campaign data
    year = campaign.x.year
    year2 = year % 100
    
    
    # add processes we are interested in
    
    process_names = [
        "data", 
        "data_mu",
        "data_tau",
        "data_e",
        #Drell-Yan
        "dy_lep",
        "dy_z2ee",
        "dy_z2mumu",
        # "dy_z2tautau",
        # "dy_z2ll",
        # "dy_lep_m10to50",
        #W + jets
        "wj",
        #diboson
        "vv", #diboson inclusive
        "ww",
        "wz",
        "zz",
        #ttbar
        "tt",#ttbar inclusive
        "tt_sl",
        "tt_dl",
        "tt_fh",
        #single top
        "st",
        #single top t-channel        
        "st_tchannel_tbar",
        "st_tchannel_t",
        # single top tW channel
        "st_twchannel_t_fh",
        "st_twchannel_t_sl",
        "st_twchannel_t_dl",
        "st_twchannel_tbar_sl",
        "st_twchannel_tbar_dl",
        "st_twchannel_tbar_fh",
    ]
    for process_name in process_names:
        # add the process
        proc = cfg.add_process(procs.get(process_name))
        #for signal datasets create special tag
        if process_name.startswith("h_"):
            proc.add_tag("signal")
            
        if proc.is_mc:
            # Updated color mapping to avoid repetition and ensure unique colors
            if proc.name == "st"                  : proc.color1 = (63, 144, 218)  # Sky Blue
            if proc.name == "st_tchannel_t"       : proc.color1 = (63, 144, 218)  # Sky Blue
            if proc.name == "st_tchannel_tbar"    : proc.color1 = (87, 144, 252)  # Dodger Blue
            if proc.name == "st_twchannel_t_sl"   : proc.color1 = (146, 218, 221) # Pale Turquoise
            if proc.name == "st_twchannel_tbar_sl": proc.color1 = (148, 164, 162) # Cadet Grey
            if proc.name == "st_twchannel_t_dl"   : proc.color1 = (169, 107, 89)  # Rosy Brown
            if proc.name == "st_twchannel_tbar_dl": proc.color1 = (200, 73, 169)  # Medium Violet Red
            if proc.name == "st_twchannel_tbar_fh": proc.color1 = (131, 45, 182)  # Amethyst
            if proc.name == "tt"                  : proc.color1 = (255, 169, 14)  # Orange
            if proc.name == "tt_sl"               : proc.color1 = (255, 169, 14)  # Orange
            if proc.name == "tt_dl"               : proc.color1 = (248, 156, 32)  # Dark Golden Rod
            if proc.name == "tt_fh"               : proc.color1 = (228, 37, 54)   # Crimson Red
            if proc.name == "vv"                  : proc.color1 = (101, 99, 100)  # Charcoal
            if proc.name == "ww"                  : proc.color1 = (101, 99, 100)  # Charcoal
            if proc.name == "zz"                  : proc.color1 = (185, 172, 112) # Olive Drab
            if proc.name == "wz"                  : proc.color1 = (122, 33, 221)  # Blue Violet
            if proc.name == "dy_lep"              : proc.color1 = (156, 156, 161) # Dark Gray
            if proc.name == "wj"                  : proc.color1 = (255, 94, 2)    # Orange Red
			
			#if proc.name == "dy_lep": proc.color1 = color=(255,204,102)
            #if proc.name == "h_ggf_tautau": proc.color1 = (51,53,204)
            #if proc.name == "wj": proc.color1 = (201,89,84)
            #if proc.name == "tt_sl": proc.color1 = (153,153,204)
            #if proc.name == "tt_dl": proc.color1 = (184,184,227)
            #if proc.name == "tt_fh": proc.color1 = (87,87,141)
            #if proc.name == "ww" : proc.color1 = (102,204,102)
            #if proc.name == "wz" : proc.color1 = (49,157,49)
            #if proc.name == "zz" : proc.color1 = (120,214,120)
            #if proc.name == "vv" : proc.color1 = (102,204,102)
            #if proc.name == "zh_htt": proc.color1 = (223,102,72)
            #if proc.name == "h_ggf_htt": proc.color1 = (51,53,204)
            #if proc.name == "vv" : proc.color1 = (102,204,102)

        # configuration of colors, labels, etc. can happen here
       

    # add datasets we need to study
    dataset_names = [
        #data
        "data_e_C",
        "data_e_D",
        "data_mu_C",
        "data_mu_D",
        "data_tau_C",
        "data_tau_D",
        #Drell-Yan
        "dy_incl",
        # "dy_lep_m10to50",
        #W+jets
        "wj_incl",
        #Diboson
        "ww",
        "wz",
        "zz",
        #ttbar
        "tt_sl",
        "tt_dl",
        "tt_fh",
        #single top t-channel
        "st_tchannel_tbar",
        "st_tchannel_t",
        # single top tW channel
        "st_twchannel_t_fh",
        "st_twchannel_t_sl",
        "st_twchannel_t_dl",
        "st_twchannel_tbar_sl",
        "st_twchannel_tbar_dl",
        "st_twchannel_tbar_fh",
        ]
    
    for dataset_name in dataset_names:
        # add the dataset
        dataset = cfg.add_dataset(campaign.get_dataset(dataset_name))
        if dataset_name.startswith("h_"):
            dataset.add_tag("signal")
        # for testing purposes, limit the number of files to 1
        for info in dataset.info.values():
            if limit_dataset_files:
                info.n_files = min(info.n_files, limit_dataset_files) #<<< REMOVE THIS FOR THE FULL DATASET

    # verify that the root process of all datasets is part of any of the registered processes
    verify_config_processes(cfg, warn=True)

  
    from httcp.config.triggers import add_triggers_run3_2022_preEE
    add_triggers_run3_2022_preEE(cfg)
    
    from httcp.config.met_filters import add_met_filters
    add_met_filters(cfg)

    # default objects, such as calibrator, selector, producer, ml model, inference model, etc
    cfg.x.default_calibrator = "main"
    cfg.x.default_selector = "main"
    cfg.x.default_producer = "main"
    cfg.x.default_weight_producer = "main"
    cfg.x.default_ml_model = None
    cfg.x.default_inference_model = "example"
    cfg.x.default_categories = ("incl",)
    cfg.x.default_variables = ("event","channel_id")

    # process groups for conveniently looping over certain processs
    # (used in wrapper_factory and during plotting)
    cfg.x.process_groups = {
        "data" : ["data_mu", "data_tau","data_e"],
        "vv"   : ["ww", "wz", "zz"],
        "tt"   : ["tt_sl","tt_dl","tt_fh"],
        "st"   : ["st_tchannel_tbar","st_tchannel_t",
               "st_twchannel_t_fh","st_twchannel_t_sl","st_twchannel_t_dl",
               "st_twchannel_tbar_sl","st_twchannel_tbar_dl","st_twchannel_tbar_fh"]
    }

    # dataset groups for conveniently looping over certain datasets
    # (used in wrapper_factory and during plotting)
    cfg.x.dataset_groups = {}

    # category groups for conveniently looping over certain categories
    # (used during plotting)
    cfg.x.category_groups = {}

    # variable groups for conveniently looping over certain variables
    # (used during plotting)
    cfg.x.variable_groups = {}

    # shift groups for conveniently looping over certain shifts
    # (used during plotting)
    cfg.x.shift_groups = {}

    # selector step groups for conveniently looping over certain steps
    # (used in cutflow tasks)
    cfg.x.selector_step_groups = {
        "default": ["json", "met_filter", "dl_res_veto", "trigger", "lepton", "jet"],
    }
    #  cfg.x.selector_step_labels = {"Initial":0, 
    #                                "Trigger": , "Muon"}
     
    # whether to validate the number of obtained LFNs in GetDatasetLFNs
    # (currently set to false because the number of files per dataset is truncated to 2)
    cfg.x.validate_dataset_lfns = False
    
    # jec configuration
    # https://twiki.cern.ch/twiki/bin/view/CMS/JECDataMC?rev=201
    jerc_postfix = ""
    if year == 2016 and campaign.x.vfp == "post":
        jerc_postfix = "APV"
    elif year == 2022 and campaign.x.tag == "post":
        jerc_postfix = "EE"
    if cfg.campaign.x.run == 2:
        jerc_campaign = f"Summer19UL{year2}{jerc_postfix}"
        jet_type = "AK4PFchs"
    elif cfg.campaign.x.run == 3:
        jerc_campaign = f"Summer{year2}{jerc_postfix}_22Sep2023"
        jet_type = "AK4PFPuppi"

    cfg.x.jec = DotDict.wrap({
        "campaign": jerc_campaign,
        "version": {2016: "V7", 2017: "V5", 2018: "V5", 2022: "V2"}[year],
        "jet_type": jet_type,
        "levels_DATA": ["L2Relative", "L2L3Residual", "L3Absolute"], # "L1L2L3Res"
        "levels_MC": ["L2Relative", "L3Absolute"], 
        "levels_for_type1_met": ["L2Relative", "L2L3Residual", "L3Absolute", "L1L2L3Res"],
        "uncertainty_sources": [
            # "AbsoluteStat",
            # "AbsoluteScale",
            # "AbsoluteSample",
            # "AbsoluteFlavMap",
            # "AbsoluteMPFBias",
            # "Fragmentation",
            # "SinglePionECAL",
            # "SinglePionHCAL",
            # "FlavorQCD",
            # "TimePtEta",
            # "RelativeJEREC1",
            # "RelativeJEREC2",
            # "RelativeJERHF",
            # "RelativePtBB",
            # "RelativePtEC1",
            # "RelativePtEC2",
            # "RelativePtHF",
            # "RelativeBal",
            # "RelativeSample",
            # "RelativeFSR",
            # "RelativeStatFSR",
            # "RelativeStatEC",
            # "RelativeStatHF",
            # "PileUpDataMC",
            # "PileUpPtRef",
            # "PileUpPtBB",
            # "PileUpPtEC1",
            # "PileUpPtEC2",
            # "PileUpPtHF",
            # "PileUpMuZero",
            # "PileUpEnvelope",
            # "SubTotalPileUp",
            # "SubTotalRelative",
            # "SubTotalPt",
            # "SubTotalScale",
            # "SubTotalAbsolute",
            # "SubTotalMC",
            "Total",
            # "TotalNoFlavor",
            # "TotalNoTime",
            # "TotalNoFlavorNoTime",
            # "FlavorZJet",
            # "FlavorPhotonJet",
            # "FlavorPureGluon",
            # "FlavorPureQuark",
            # "FlavorPureCharm",
            # "FlavorPureBottom",
            # # "TimeRunA",
            # # "TimeRunB",
            # # "TimeRunC",
            # # "TimeRunD",
            # "CorrelationGroupMPFInSitu",
            # "CorrelationGroupIntercalibration",
            # "CorrelationGroupbJES",
            # "CorrelationGroupFlavor",
            # "CorrelationGroupUncorrelated",
        ],
    })
    
    ###############################################################################################
    # met settings
    ################################################################################################

    # name of the MET phi correction set
    # (used in the met_phi calibrator)
    #cfg.x.met_phi_correction_set = r"{variable}_metphicorr_pfmet_{data_source}"
    
    ###############################################################################################
    # JER settings
    ################################################################################################
    
    # # JER
    # # https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution?rev=107
    # # TODO: get jerc working for Run3
    # cfg.x.jer = DotDict.wrap({
    #     "campaign": jerc_campaign,
    #     "version": {2016: "JRV3", 2017: "JRV2", 2018: "JRV2", 2022: "JRV1"}[year],
    #     "jet_type": jet_type,
    # })

    # # JEC uncertainty sources propagated to btag scale factors
    # # (names derived from contents in BTV correctionlib file)
    # cfg.x.btag_sf_jec_sources = [
    #     "",  # total
    #     "Absolute",
    #     "AbsoluteMPFBias",
    #     "AbsoluteScale",
    #     "AbsoluteStat",
    #     f"Absolute_{year}",
    #     "BBEC1",
    #     f"BBEC1_{year}",
    #     "EC2",
    #     f"EC2_{year}",
    #     "FlavorQCD",
    #     "Fragmentation",
    #     "HF",
    #     f"HF_{year}",
    #     "PileUpDataMC",
    #     "PileUpPtBB",
    #     "PileUpPtEC1",
    #     "PileUpPtEC2",
    #     "PileUpPtHF",
    #     "PileUpPtRef",
    #     "RelativeBal",
    #     "RelativeFSR",
    #     "RelativeJEREC1",
    #     "RelativeJEREC2",
    #     "RelativeJERHF",
    #     "RelativePtBB",
    #     "RelativePtEC1",
    #     "RelativePtEC2",
    #     "RelativePtHF",
    #     "RelativeSample",
    #     f"RelativeSample_{year}",
    #     "RelativeStatEC",
    #     "RelativeStatFSR",
    #     "RelativeStatHF",
    #     "SinglePionECAL",
    #     "SinglePionHCAL",
    #     "TimePtEta",
    # ]

    # lumi values in inverse pb
    # https://twiki.cern.ch/twiki/bin/view/CMS/LumiRecommendationsRun2?rev=2#Combination_and_correlations
    #https://twiki.cern.ch/twiki/bin/viewauth/CMS/PdmVRun3Analysis#DATA_AN2
    #Only F and G eras
    cfg.x.luminosity = Number(7980, {
        "lumi_13p6TeV_2022": 0.014j,
        
    })
    

    # names of muon correction sets and working points
    # (used in the muon producer)
   
  
    cfg.x.deep_tau = DotDict.wrap({
        "tagger": "DeepTau2018v2p5",
        "vs_e"          : "VVLoose",
        "vs_mu"         : "Tight",
        "vs_jet"        : "Medium",
        "vs_e_jet_wps"  : {'VVVLoose'   : 1,
                           'VVLoose'    : 2,
                           'VLoose'     : 3,
                           'Loose'      : 4,
                           'Medium'     : 5,
                           'Tight'      : 6,
                           'VTight'     : 7,
                           'VVTight'    : 8},
        "vs_mu_wps"     : {'VLoose' : 1,
                           'Loose'  : 2,
                           'Medium' : 3,
                           'Tight'  : 4}
        })

    cfg.x.btag_working_points = DotDict.wrap(
        {
            2016 : {
                "deepjet": { #TODO: make a link to this numbers
                    "loose": 0.0532,
                    "medium": 0.3040,
                    "tight": 0.7476,
                },
                "deepcsv": {
                    "loose": 0.1355,
                    "medium": 0.4506,
                    "tight": 0.7738,
                },
            },
            2022 : {
                "deepjet" : { #https://btv-wiki.docs.cern.ch/ScaleFactors/Run3Summer22/
                    "loose": 0.0583,
                    "medium": 0.3086,
                    "tight": 0.7183,
                }
            }
                    
                
        },
    )
    
    def tag_caster(campaign: od.Campaign) -> str:
        #Helper function to cast campaign tags to the tags used in POG groups for the scale factors
        year = campaign.x.year
        tag = campaign.x.tag
        out_tag = ''
        if year in [2017,2018]  : out_tag = 'UL'
        elif tag == "preEE"     : out_tag = "Summer22"
        elif tag == "postEE"    : out_tag = "Summer22EE"
        elif tag == "preBpix"   : out_tag = "Summer23"
        elif tag == "postBpix"  : out_tag = "Summer23BPix"
        elif tag == "preVFP"    : out_tag = "preVFP_UL"
        elif tag == "postVFP"   : out_tag = "postVFP_UL"
        return out_tag

    import os
    tag = tag_caster(campaign)
    ### Jet energy correction and resolution configuration ###
    cfg.x.jec = DotDict.wrap({
        "campaign": f"{tag}_22Sep2023",
        "version":  "V2",
        "jet_type": "AK4PFPuppi",
        "levels": ["L1FastJet", "L2Relative", "L2L3Residual", "L3Absolute"],
        "levels_for_type1_met": ["L1FastJet"],
        "uncertainty_sources": [
            # "AbsoluteStat",
            # "AbsoluteScale",
            # "AbsoluteSample",
            # "AbsoluteFlavMap",
            # "AbsoluteMPFBias",
            # "Fragmentation",
            # "SinglePionECAL",
            # "SinglePionHCAL",
            # "FlavorQCD",
            # "TimePtEta",
            # "RelativeJEREC1",
            # "RelativeJEREC2",
            # "RelativeJERHF",
            # "RelativePtBB",
            # "RelativePtEC1",
            # "RelativePtEC2",
            # "RelativePtHF",
            # "RelativeBal",
            # "RelativeSample",
            # "RelativeFSR",
            # "RelativeStatFSR",
            # "RelativeStatEC",
            # "RelativeStatHF",
            # "PileUpDataMC",
            # "PileUpPtRef",
            # "PileUpPtBB",
            # "PileUpPtEC1",
            # "PileUpPtEC2",
            # "PileUpPtHF",
            # "PileUpMuZero",
            # "PileUpEnvelope",
            # "SubTotalPileUp",
            # "SubTotalRelative",
            # "SubTotalPt",
            # "SubTotalScale",
            # "SubTotalAbsolute",
            # "SubTotalMC",
            "Total",
            # "TotalNoFlavor",
            # "TotalNoTime",
            # "TotalNoFlavorNoTime",
            # "FlavorZJet",
            # "FlavorPhotonJet",
            # "FlavorPureGluon",
            # "FlavorPureQuark",
            # "FlavorPureCharm",
            # "FlavorPureBottom",
            "CorrelationGroupMPFInSitu",
            "CorrelationGroupIntercalibration",
            "CorrelationGroupbJES",
            "CorrelationGroupFlavor",
            "CorrelationGroupUncorrelated",
        ],
    })
    
    cfg.x.jer = DotDict.wrap({
        "campaign": f"{tag}_22Sep2023",
        "version": "JRV1",
        "jet_type": "AK4PFPuppi",
        })
    
    corr_dir = os.path.join(os.environ.get('CF_REPO_BASE'), "httcp/corrections/")
    jsonpog_dir = os.path.join(os.environ.get('CF_REPO_BASE'), "modules/jsonpog-integration/POG/")


    cfg.x.external_files = DotDict.wrap({
        # lumi files
        "lumi": {
            "golden": (f"{corr_dir}Cert_Collisions2022_355100_362760_Golden.json", "v1"),  # https://cms-service-dqmdc.web.cern.ch/CAF/certification/Collisions18/13TeV/Legacy_2018/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt
            "normtag": (f"{corr_dir}normtag_PHYSICS.json", "v1"), #/cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags
        },
        "pu_sf": (f"{jsonpog_dir}LUM/{cfg.x.year}{tag}/puWeights.json.gz", "v1"),
        "muon_correction" : f"{jsonpog_dir}MUO/{cfg.x.year}{tag}/muon_Z.json.gz",
        #"electron_correction" : (f"{jsonpog_dir}EGM/{cfg.x.year}{tag}/electron.json.gz", "v2")
        #"tau_correction"  : f"{jsonpog_dir}TAU/{cfg.x.year}{tag}/tau.json.gz", 
        "tau_correction"  : f"{corr_dir}tau_DeepTau2018v2p5_2022_preEE.json.gz", #FIXME: this sf json is not from the jsonpog-integration dir!
        "jet_jerc"  : (f"{jsonpog_dir}JME/{cfg.x.year}{tag}/jet_jerc.json.gz", "v2"),
        "jet_veto_map"  : (f"{jsonpog_dir}JME/{cfg.x.year}{tag}/jetvetomaps.json.gz", "v2"),
        #"met_phi_corr": (f"{jsonpog_dir}JME/{cfg.x.year}{tag}/met{cfg.x.year}.json.gz", "v2"), #FIXME: there is no json present in the jsonpog-integration for this year, I retrieve the json frm: https://cms-talk.web.cern.ch/t/2022-met-xy-corrections/53414/2 but it seems corrupted
    })

    # target file size after MergeReducedEvents in MB
    cfg.x.reduced_file_size = 512.0
    
    from httcp.config.variables import keep_columns
    keep_columns(cfg)
 
    cfg.add_shift(name="nominal", id=0)

    cfg.add_shift(name="tau_up", id=1, type="shape")
    cfg.add_shift(name="tau_down", id=2, type="shape")
    add_shift_aliases(cfg, "tau", {"tau_weight": "tau_weight_{direction}"})
    
    cfg.add_shift(name="mu_up", id=3, type="shape")
    cfg.add_shift(name="mu_down", id=4, type="shape")
    add_shift_aliases(cfg, "mu", {"muon_weight": "muon_weight_{direction}"})
    
    cfg.add_shift(name="ts_up", id=5, type="shape") #cp-even
    cfg.add_shift(name="ts_down", id=7, type="shape") #cp-odd
    add_shift_aliases(cfg, "ts", {"tauspinner_weight": "tauspinner_weight_{direction}"})
    # event weight columns as keys in an OrderedDict, mapped to shift instances they depend on
    get_shifts = functools.partial(get_shifts_from_sources, cfg)   

    # thisdir = os.path.dirname(os.path.abspath(__file__))
    
    # with open(os.path.join(thisdir, "jec_sources.yaml"), "r") as f:
    #     all_jec_sources = yaml.load(f, yaml.Loader)["names"]

    # for jec_source in cfg.x.jec["uncertainty_sources"]:
    #     idx = all_jec_sources.index(jec_source)
    #     cfg.add_shift(
    #         name=f"jec_{jec_source}_up",
    #         id=5000 + 2 * idx,
    #         type="shape",
    #         tags={"jec"},
    #         aux={"jec_source": jec_source},
    #     )
    #     cfg.add_shift(
    #         name=f"jec_{jec_source}_down",
    #         id=5001 + 2 * idx,
    #         type="shape",
    #         tags={"jec"},
    #         aux={"jec_source": jec_source},
    #     )
    #     add_shift_aliases(
    #         cfg,
    #         f"jec_{jec_source}",
    #         {"Jet.pt": "Jet.pt_{name}", "Jet.mass": "Jet.mass_{name}"},
    #     )

    #     if jec_source in ["Total", *cfg.x.btag_sf_jec_sources]:
    #         # when jec_source is a known btag SF source, add aliases for btag weight column
    #         add_shift_aliases(
    #             cfg,
    #             f"jec_{jec_source}",
    #             {
    #                 "btag_weight": f"btag_weight_jec_{jec_source}_" + "{direction}",
    #                 "normalized_btag_weight": f"normalized_btag_weight_jec_{jec_source}_" + "{direction}",
    #                 "normalized_njet_btag_weight": f"normalized_njet_btag_weight_jec_{jec_source}_" + "{direction}",
    #             },
    #         )

    # cfg.add_shift(name="jer_up", id=6000, type="shape", tags={"jer"})
    # cfg.add_shift(name="jer_down", id=6001, type="shape", tags={"jer"})
    # add_shift_aliases(cfg, "jer", {"Jet.pt": "Jet.pt_{name}", "Jet.mass": "Jet.mass_{name}"})
  
    

    
    # versions per task family, either referring to strings or to callables receving the invoking
    # task instance and parameters to be passed to the task family
    def set_version(cls, inst, params):
        # per default, use the version set on the command line
        version = inst.version 
        return version if version else 'dev'
            
        
    cfg.x.versions = {
        "cf.CalibrateEvents"    : set_version,
        "cf.SelectEvents"       : set_version,
        "cf.MergeSelectionStats": set_version,
        "cf.MergeSelectionMasks": set_version,
        "cf.ReduceEvents"       : set_version,
        "cf.MergeReductionStats": set_version,
        "cf.MergeReducedEvents" : set_version,
    }
    # channels
    cfg.add_channel(name="etau",   id=1)
    cfg.add_channel(name="mutau",  id=2)
    #cfg.add_channel(name="emu"  ,  id=3)
    cfg.add_channel(name="tautau", id=4)
    
    cfg.x.ch_objects = DotDict.wrap({
        'etau'   : {'lep0' : 'Electron',
                    'lep1' : 'Tau'},
        'mutau'  : {'lep0' : 'Muon',
                    'lep1' : 'Tau'},
        'tautau' : {'lep0' : 'Tau',
                    'lep1' : 'Tau'},
    })
    
    if cfg.campaign.x("custom").get("creator") == "desy":  
        def get_dataset_lfns(dataset_inst: od.Dataset, shift_inst: od.Shift, dataset_key: str) -> list[str]:
            # destructure dataset_key into parts and create the lfn base directory
            print(f"Creating custom get_dataset_lfns for {config_name}")   
            try:
               basepath = cfg.campaign.x("custom").get("location")
            except:
                print("Did not find any basebath in the campaigns")
                basepath = "" 
            lfn_base = law.wlcg.WLCGDirectoryTarget(
                f"{basepath}{dataset_key}",
                fs="wlcg_fs_eos",
            )
            print(f"lfn basedir:{lfn_base}")
            # loop though files and interpret paths as lfns
            return [
                lfn_base.child(basename, type="f").path
                for basename in lfn_base.listdir(pattern="*.root")
            ]
        # define the lfn retrieval function
        cfg.x.get_dataset_lfns = get_dataset_lfns
        # define a custom sandbox
        cfg.x.get_dataset_lfns_sandbox = dev_sandbox("bash::$CF_BASE/sandboxes/cf.sh")
        # define custom remote fs's to look at
        cfg.x.get_dataset_lfns_remote_fs =  lambda dataset_inst: "wlcg_fs_eos"
        
    # add categories using the "add_category" tool which adds auto-generated ids
    from httcp.config.categories import add_categories
    add_categories(cfg)
        
    from httcp.config.variables import add_variables
    add_variables(cfg)
    
    
