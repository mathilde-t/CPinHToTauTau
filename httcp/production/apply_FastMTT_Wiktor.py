# coding: utf-8

"""
Column production methods related to higher-level features.
"""
import functools

import law
import order as od
from typing import Optional
from columnflow.production import Producer, producer
from columnflow.production.util import attach_coffea_behavior

from columnflow.util import maybe_import
from columnflow.columnar_util import EMPTY_FLOAT, Route, set_ak_column, remove_ak_column
from columnflow.columnar_util import optional_column as optional

from modules.ClassicsSVfit.python.FastMTT import FastMTT

np = maybe_import("numpy")
ak = maybe_import("awkward")
pd = maybe_import("pandas")
coffea = maybe_import("coffea")
maybe_import("coffea.nanoevents.methods.nanoaod")

# helpers
set_ak_column_f32 = functools.partial(set_ak_column, value_type=np.float32)
set_ak_column_i32 = functools.partial(set_ak_column, value_type=np.int32)

logger = law.logger.get_logger(__name__)


@producer(
    uses={
        # nano columns
        "hcand.*", "channel_id",
        "PuppiMET.pt", "PuppiMET.phi", "PuppiMET.covXX", "PuppiMET.covXY", "PuppiMET.covYY",
        optional("GenTau.pt"), optional("GenTau.eta"), optional("GenTau.phi"), optional("GenTau.mass"),
    },
    produces={
        # new columns
        "hcand.pt_fastMTT_W", "hcand.eta_fastMTT_W", "hcand.phi_fastMTT_W", "hcand.mass_fastMTT_W",
        "hcand_invm_fastMTT_W",
        "hcand.px_fastMTT", "hcand.py_fastMTT", "hcand.pz_fastMTT",

        ### gen_reco resolutions
        "hcand.px_fastMTT_W_resolution_gen_reco","hcand.py_fastMTT_W_resolution_gen_reco","hcand.pz_fastMTT_W_resolution_gen_reco",
        "hcand.pt_fastMTT_W_resolution_gen_reco_relative","hcand.pt_fastMTT_W_resolution_gen_reco_absolute",
        "hcand.eta_fastMTT_W_resolution_gen_reco","hcand.phi_fastMTT_W_resolution_gen_rec",
        ### gen_fastmtt resolutions
        "hcand.px_fastMTT_W_resolution_gen_fastmtt","hcand.py_fastMTT_W_resolution_gen_fastmtt","hcand.pz_fastMTT_W_resolution_gen_fastmtt",
        "hcand.pt_fastMTT_W_resolution_gen_fastmtt_relative","hcand.pt_fastMTT_W_resolution_gen_fastmtt_absolute",
        "hcand.eta_fastMTT_W_resolution_gen_fastmtt","hcand.phi_fastMTT_W_resolution_gen_fastmt",
        ### reco_fastmtt resolutions
        "hcand.px_fastMTT_W_resolution_reco_fastmtt","hcand.py_fastMTT_W_resolution_reco_fastmtt","hcand.pz_fastMTT_W_resolution_reco_fastmtt",
        "hcand.pt_fastMTT_W_resolution_reco_fastmtt_relative","hcand.pt_fastMTT_W_resolution_reco_fastmtt_absolute",
        "hcand.eta_fastMTT_W_resolution_reco_fastmtt","hcand.phi_fastMTT_W_resolution_reco_fastmt",

    },
)
def apply_fastMTT_Wiktor(
        self: Producer, 
        events: ak.Array,
        **kwargs
) -> ak.Array:

    events = ak.Array(events, behavior=coffea.nanoevents.methods.nanoaod.behavior)
    hcand_ = ak.with_name(events.hcand, "PtEtaPhiMLorentzVector")
    hmass  = (hcand_[:,:1] + hcand_[:,1:2]).mass
    
    etau_id   = self.config_inst.get_channel("etau").id
    mutau_id  = self.config_inst.get_channel("mutau").id
    tautau_id = self.config_inst.get_channel("tautau").id
    
    # decay_type:
    #1 - TauToHad
    #2 - TauToElec
    #3 - TauToMu
    type_1 = ak.ones_like(hmass)
    type_2 = 2 * type_1
    type_3 = 3 * type_1
    
    pt1    = ak.to_numpy(events.hcand.pt[:,0:1])
    eta1   = ak.to_numpy(events.hcand.eta[:,0:1])
    phi1   = ak.to_numpy(events.hcand.phi[:,0:1])
    mass1  = ak.to_numpy(events.hcand.mass[:,0:1])
    dm1    = events.hcand.decayMode[:,0:1]
    dm1_dummy = ak.values_astype(-1 * ak.ones_like(dm1), np.int32)
    dm1    = ak.to_numpy(ak.where(events.channel_id < 4, dm1_dummy, dm1))
    type1  = ak.to_numpy(ak.where(events.channel_id == etau_id,
                                  type_2,
                                  ak.where(events.channel_id == mutau_id,
                                           type_3,
                                           type_1)
                                  )
                         )
    
    pt2    = ak.to_numpy(events.hcand.pt[:,1:2])
    eta2   = ak.to_numpy(events.hcand.eta[:,1:2])
    phi2   = ak.to_numpy(events.hcand.phi[:,1:2])
    mass2  = ak.to_numpy(events.hcand.mass[:,1:2])
    dm2    = ak.to_numpy(events.hcand.decayMode[:,1:2])
    type2  = ak.to_numpy(type_1)
    
    metpt  = ak.to_numpy(events.PuppiMET.pt[:,None])
    metphi = ak.to_numpy(events.PuppiMET.phi[:,None])
    metcovxx = ak.to_numpy(events.PuppiMET.covXX[:,None])
    metcovxy = ak.to_numpy(events.PuppiMET.covXY[:,None])
    metcovyy = ak.to_numpy(events.PuppiMET.covYY[:,None])

    higgs_mass = ak.to_numpy(hmass)

    gen_pt1    = ak.to_numpy(events.GenTau.pt[:,0:1])
    gen_eta1   = ak.to_numpy(events.GenTau.eta[:,0:1])
    gen_phi1   = ak.to_numpy(events.GenTau.phi[:,0:1])
    gen_mass1  = ak.to_numpy(events.GenTau.mass[:,0:1])
    gen_pt2    = ak.to_numpy(events.GenTau.pt[:,1:2])
    gen_eta2   = ak.to_numpy(events.GenTau.eta[:,1:2])
    gen_phi2   = ak.to_numpy(events.GenTau.phi[:,1:2])
    gen_mass2  = ak.to_numpy(events.GenTau.mass[:,1:2])

    # prepare inputs properly
    metcov = np.concatenate((metcovxx,metcovxy,metcovxy,metcovyy), axis=1)
    metcov = np.reshape(metcov, (metcov.shape[0], 2, 2))
    measuredTauLeptons = np.concatenate((dm1,pt1,eta1,phi1,mass1,type1,
                                         dm2,pt2,eta2,phi2,mass2,type2), 
                                         axis=1)
    measuredTauLeptons = np.reshape(measuredTauLeptons, (measuredTauLeptons.shape[0], 2, 6))
    metx   = (metpt * np.cos(metphi)).reshape(-1)
    mety   = (metpt * np.sin(metphi)).reshape(-1)
    
    
    # Launch FastMTT
    fMTT = FastMTT()
    #You can choose to plot likelihood for one of the events. -1 means no plot.
    fMTT.WhichLikelihoodPlot = -1
    #You can also choose to calculate uncertainties by:
    fMTT.CalculateUncertainties = True

    fMTT.run(measuredTauLeptons, metx, mety, metcov)

    mass_h = fMTT.mass
    mass_h = ak.from_regular(ak.Array(mass_h.reshape(mass_h.shape[0],1)))
    
    p4_h1   = fMTT.tau1P4
    px_h1   = ak.from_regular(ak.Array(p4_h1[:,0:1]))
    py_h1   = ak.from_regular(ak.Array(p4_h1[:,1:2]))
    pz_h1   = ak.from_regular(ak.Array(p4_h1[:,2:3]))
    energy_h1 = ak.from_regular(ak.Array(p4_h1[:,3:4]))
    p4_h1_reg = ak.zip(
        {
            "x": px_h1, "y": py_h1, "z": pz_h1, "t": energy_h1,
        },
        with_name="LorentzVector",
        behavior=coffea.nanoevents.methods.vector.behavior,
    )
    
    p4_h2  = fMTT.tau2P4
    px_h2   = ak.from_regular(ak.Array(p4_h2[:,0:1]))
    py_h2   = ak.from_regular(ak.Array(p4_h2[:,1:2]))
    pz_h2   = ak.from_regular(ak.Array(p4_h2[:,2:3]))
    energy_h2 = ak.from_regular(ak.Array(p4_h2[:,3:4]))
    p4_h2_reg = ak.zip(
        {
            "x": px_h2, "y": py_h2, "z": pz_h2, "t": energy_h2,
        },
        with_name="LorentzVector",
        behavior=coffea.nanoevents.methods.vector.behavior,
    )

    hcand_px_fastMTT = ak.concatenate([px_h1, px_h2], axis=1)
    hcand_py_fastMTT = ak.concatenate([py_h1, py_h2], axis=1)
    hcand_pz_fastMTT = ak.concatenate([pz_h1, pz_h2], axis=1)
    hcand_pt_fastMTT  = ak.concatenate([p4_h1_reg.pt, p4_h2_reg.pt], axis=1)
    hcand_eta_fastMTT = ak.concatenate([p4_h1_reg.eta, p4_h2_reg.eta], axis=1)
    hcand_phi_fastMTT = ak.concatenate([p4_h1_reg.phi, p4_h2_reg.phi], axis=1)
    h1_mass = p4_h1_reg.mass
    logger.critical("check regressed mass for leptonic channel, many nan values found for etau channel")
    h1_mass = ak.where(events.channel_id == etau_id,
                       ak.nan_to_num(h1_mass, nan=0.00051),
                       ak.where(events.channel_id == mutau_id,
                                ak.nan_to_num(h1_mass, nan=0.10566),
                                h1_mass))

    hcand_mass_fastMTT = ak.concatenate([h1_mass, p4_h2_reg.mass], axis=1)

    ## only plot the first 300 events to be coherent with IC's
    ## for test purposes only !
    #mask = np.arange(len(events)) < 300 # initialise Mask : True für die ersten 300 Events, False für den Rest
    #mask = ak.Array(mask)
    #zero_array = lambda x: ak.zeros_like(x)
    #
    #hcand_pt_fastMTT = ak.where(mask, hcand_pt_fastMTT, zero_array(hcand_pt_fastMTT))
    #hcand_eta_fastMTT = ak.where(mask, hcand_eta_fastMTT, zero_array(hcand_eta_fastMTT))
    #hcand_phi_fastMTT = ak.where(mask, hcand_phi_fastMTT, zero_array(hcand_phi_fastMTT))
    #hcand_mass_fastMTT = ak.where(mask, hcand_mass_fastMTT, zero_array(hcand_mass_fastMTT))
    #mass_h = ak.where(mask, mass_h, zero_array(mass_h))

    events = set_ak_column(events, "hcand.px_fastMTT",   hcand_px_fastMTT)
    events = set_ak_column(events, "hcand.py_fastMTT",   hcand_py_fastMTT)
    events = set_ak_column(events, "hcand.pz_fastMTT",   hcand_pz_fastMTT)
    events = set_ak_column(events, "hcand.pt_fastMTT_W",   hcand_pt_fastMTT)
    events = set_ak_column(events, "hcand.eta_fastMTT_W",  hcand_eta_fastMTT)
    events = set_ak_column(events, "hcand.phi_fastMTT_W",  hcand_phi_fastMTT)
    events = set_ak_column(events, "hcand.mass_fastMTT_W", hcand_mass_fastMTT)

    events = set_ak_column_f32(events, "hcand_invm_fastMTT_W", mass_h)



    #################################
    ####       resolution         ###
    #################################

    def compute_resolutions(events):
        hcands = [("hcand1", "h1", 1), ("hcand2", "h2", 2)]
        variables = ["px", "py", "pz", "pt", "eta", "phi"]
        comparisons = {
            "gen_reco": ("gen_{var}{num}", "{var}{num}"),
            "gen_fastmtt": ("gen_{var}{num}", "p4_{h_var}_reg.{var}"),
            "reco_fastmtt": ("p4_{h_var}_reg.{var}", "{var}{num}")
        }
    
        resolutions = {}
        for hcand, h_var, num in hcands:
            for var in variables:
                for comp, (num_fmt, denom_fmt) in comparisons.items():
                    num_key = num_fmt.format(var=var, num=num, h_var=h_var)
                    denom_key = denom_fmt.format(var=var, num=num, h_var=h_var)
                    
                    num_val = events[num_key] if num_key in events.fields else None
                    denom_val = events[denom_key] if denom_key in events.fields else None
                    
                    if num_val is not None and denom_val is not None:
                        key = f"resolution_{hcand}_{var}_fastMTT_resolution_{comp}"
                        #resolutions[key] = (denom_val - num_val) / denom_val if var == "pt" else denom_val - num_val
                        if "pt" in var:
                            resolutions[key] = {
                                "relative": (denom_val - num_val) / denom_val,
                                "absolute": denom_val - num_val
                            }
                        else:
                            resolutions[key] = denom_val - num_val

    
    ## Concatenate results
    #for comp in comparisons.keys():
    #    for var in variables:
    #        key_list = [resolutions[f"resolution_{hcand}_{var}_fastMTT_resolution_{comp}"] for hcand, _, _ in hcands]
    #        concatenated = ak.concatenate(key_list, axis=1)
    #        events = set_ak_column(events, f"hcand.{var}_fastMTT_W_resolution_{comp}", concatenated)

    for comp in comparisons.keys():
        for var in variables:
            if "pt" in var:
                key_list_rel = [resolutions[f"resolution_{hcand}_{var}_fastMTT_resolution_{comp}"]["relative"] for hcand, _, _ in hcands]
                key_list_abs = [resolutions[f"resolution_{hcand}_{var}_fastMTT_resolution_{comp}"]["absolute"] for hcand, _, _ in hcands]

                concatenated_rel = ak.concatenate(key_list_rel, axis=1)
                concatenated_abs = ak.concatenate(key_list_abs, axis=1)

                events = set_ak_column(events, f"hcand.{var}_fastMTT_W_resolution_{comp}_relative", concatenated_rel)
                events = set_ak_column(events, f"hcand.{var}_fastMTT_W_resolution_{comp}_absolute", concatenated_abs)
            else:
                key_list = [resolutions[f"resolution_{hcand}_{var}_fastMTT_resolution_{comp}"] for hcand, _, _ in hcands]
                concatenated = ak.concatenate(key_list, axis=1)
                events = set_ak_column(events, f"hcand.{var}_fastMTT_W_resolution_{comp}", concatenated)



    return events