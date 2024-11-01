# coding: utf-8

"""
Collection of helpers
"""

from __future__ import annotations


import law
import order as od
from typing import Any
from columnflow.util import maybe_import
from columnflow.columnar_util import ArrayFunction, deferred_column
from columnflow.util import DotDict

np = maybe_import("numpy")
ak = maybe_import("awkward")
coffea = maybe_import("coffea")
maybe_import("coffea.nanoevents.methods.nanoaod")


@deferred_column
def IF_NANO_V9(self, func: ArrayFunction) -> Any | set[Any]:
    return self.get() if func.config_inst.campaign.x.version == 9 else None


@deferred_column
def IF_NANO_V11(self, func: ArrayFunction) -> Any | set[Any]:
    return self.get() if func.config_inst.campaign.x.version >= 10 else None

# lambda function to get 4-vector of a lepton
def get_lep_p4(part): return ak.zip({f"{var}": part[var] for var in ['pt', 'eta', 'phi', 'mass']},
                                    with_name="PtEtaPhiMLorentzVector",
                                    behavior=coffea.nanoevents.methods.vector.behavior)

# lambda function to get 4-vector of impact parameter from the particle objects
# Zeroth component of IP vector is set to zero by definition that can be found here: https://www.mdpi.com/2218-1997/8/5/256
def get_ip_p4(part): return ak.zip({f'{var}': part[f'IP{var}']for var in ['x', 'y', 'z']} | {'t': ak.zeros_like(part.IPx)},
                                   with_name="LorentzVector",
                                   behavior=coffea.nanoevents.methods.vector.behavior)# # lambda function to get 4-vector from the particle objects


def transverse_mass(lepton: ak.Array, met: ak.Array) -> ak.Array:
    dphi_lep_met = lepton.delta_phi(met)
    mt = np.sqrt(2 * lepton.pt * met.pt * (1 - np.cos(dphi_lep_met)))
    return mt


def trigger_object_matching(
    vectors1: ak.Array,
    vectors2: ak.Array,
    threshold: float = 0.5,
    axis: int = 2,
) -> ak.Array:
    """
    Helper to check per object in *vectors1* if there is at least one object in *vectors2* that
    leads to a delta R metric below *threshold*. The final reduction is applied over *axis* of the
    resulting metric table containing the full combinatorics. When *return_all_matches* is *True*,
    the matrix with all matching decisions is returned as well.
    """
    # delta_r for all combinations
    dr = vectors1.metric_table(vectors2)
    # check per element in vectors1 if there is at least one matching element in vectors2
    any_match = (dr < threshold)
    return any_match


def hlt_path_fired(events,dictionary):
    HLT_path_fired = -1*ak.ones_like(ak.local_index(events.event), dtype=np.int64)
    if len(dictionary) > 0:
        # Start by collecting all the arrays from the dictionary
        # Collect arrays from the dictionary
        array_list = []
        for key in dictionary.keys():
            array_list.append(dictionary[key])

        # Stack all the arrays along the second axis (axis=1)
        hlt_condition_values_concat = np.stack(array_list, axis=1)
        HLT_path_fired = ak.fill_none(ak.max(hlt_condition_values_concat, axis=-1),-1)
    return HLT_path_fired 
    
# Define the fill_missing function
def fill_missing(arr, fill_value):
    return ak.fill_none(arr, fill_value)
   
def get_dataset_lfns(
        dataset_inst: od.Dataset,
        shift_inst: od.Shift,
        dataset_key: str,
) -> list[str]:
    # destructure dataset_key into parts and create the lfn base directory
    lfn_base = law.wlcg.WLCGDirectoryTarget(
        dataset_key,
        fs=f"local",
    )
    # loop though files and interpret paths as lfns
    paths = [lfn_base.child(basename, type="f").path for basename in lfn_base.listdir(pattern="*.root")]

    return paths


def getGenTauDecayMode(prod: ak.Array):
    pids = prod.pdgId

    is_ele  = np.abs(pids) == 11
    is_muon = np.abs(pids) == 13
    is_charged = ((np.abs(pids) == 211) | (np.abs(pids) == 321))
    is_neutral = ((pids == 111) | (pids == 311) | (pids == 130) | (pids == 310))

    edecay = ak.sum(is_ele,  axis=-1) > 0
    mdecay = ak.sum(is_muon, axis=-1) > 0
    hdecay = (ak.sum(is_charged, axis=-1) > 0) | (ak.sum(is_neutral, axis=-1) >= 0)

    Nc = ak.sum(is_charged, axis=-1)
    Np = ak.sum(is_neutral, axis=-1)

    dm = ak.where(edecay, 
                  -1, 
                  ak.where(mdecay, 
                           -2, 
                           ak.where(hdecay, 
                                    (5 * (Nc - 1) + Np),
                                    -9)
                       )
              )

    return dm


def has_pt_greater_equal(obj_1, obj_2, pt_offset):
    """
    Checks if the pt of any object in vector1 is greater than or equal to 
    the pt of any object in vector2 plus an offset.
    
    Parameters:
        vector1 (awkward.Array): The pt values of the first vector (nested).
        vector2 (awkward.Array): The pt values of the second vector (nested).
        offset (float): The value to add to the pt of objects in vector2.
    
    Returns:
        awkward.Array: A boolean array with True if at least one object in vector1 
                       has pt greater than or equal to pt in vector2 + offset, otherwise False.
    """
    # Cartesian product to get all pairs
    pairs = ak.cartesian([obj_1, obj_2], axis=-1, nested=True)
    obj_1_br, obj_2_br = ak.unzip(pairs)
    return obj_1_br.pt >= (obj_2_br.pt + pt_offset)

def hlt_path_matching(events, triggers, pair_objects):
    # Initialize masks and dictionaries
    false_mask = ak.zeros_like(ak.local_index(events.event), dtype=np.bool_)
    single_electron_triggered = false_mask
    cross_electron_triggered  = false_mask
    single_muon_triggered = false_mask
    cross_muon_triggered  = false_mask
    cross_tau_triggered  = false_mask
    
    hlt_path_fired_e     = {}
    hlt_path_fired_etau  = {}
    hlt_path_fired_mu    = {}
    hlt_path_fired_mutau = {}
    hlt_path_fired_tau   = {}

    # Perform each lepton election step separately per trigger
    for trigger, trigger_fired, leg_masks in triggers.x.trigger_data:
        
        is_single_el = trigger.has_tag("single_e")
        is_cross_el  = trigger.has_tag("cross_e_tau")
        is_single_mu = trigger.has_tag("single_mu")
        is_cross_mu  = trigger.has_tag("cross_mu_tau")
        is_cross_tau = trigger.has_tag("cross_tau_tau")
        
        if is_single_mu or is_cross_mu:
            muons = pair_objects.mutau.lep0
            taus = pair_objects.mutau.lep1
            
            if is_single_mu:
                assert trigger.n_legs == len(leg_masks) == 1
                assert abs(trigger.legs[0].pdg_id) == 13
                dr_matching = trigger_object_matching(muons, events.TrigObj[leg_masks[0]])
                single_mu_matches_leg0 = dr_matching
                single_mu_matches_leg0 = ak.any(ak.flatten(single_mu_matches_leg0, axis=-1), axis=1)
                single_muon_triggered = ak.where(trigger_fired & is_single_mu, True, single_muon_triggered)
                hlt_path_fired_mu[trigger.hlt_field] = ak.where(single_mu_matches_leg0, trigger.id, -1)
            
            elif is_cross_mu:
                assert trigger.n_legs == len(leg_masks) == 2
                assert abs(trigger.legs[0].pdg_id) == 13
                assert abs(trigger.legs[1].pdg_id) == 15
                
                dr_matching_mu = trigger_object_matching(muons, events.TrigObj[leg_masks[0]])
                cross_mu_matches_leg0 = dr_matching_mu
                dr_matching_tau = trigger_object_matching(taus, events.TrigObj[leg_masks[1]])
                cross_mu_tau_matches_leg1 = dr_matching_tau
                cross_mu_tau_matched = (ak.any(ak.flatten(cross_mu_matches_leg0, axis=-1), axis=1) & 
                                        ak.any(ak.flatten(cross_mu_tau_matches_leg1, axis=-1), axis=1))
                cross_muon_triggered = ak.where(trigger_fired & is_cross_mu, True, cross_muon_triggered)
                hlt_path_fired_mutau[trigger.hlt_field] = ak.where(cross_mu_tau_matched, trigger.id, -1)

        if is_single_el or is_cross_el:
            electrons =  pair_objects.etau.lep0
            taus =  pair_objects.etau.lep1

            if is_single_el:
                assert trigger.n_legs == len(leg_masks) == 1
                assert abs(trigger.legs[0].pdg_id) == 11
                dr_matching_e = trigger_object_matching(electrons, events.TrigObj[leg_masks[0]])
                single_e_matches_leg0 =  dr_matching_e
                single_e_matches_leg0 = ak.any(ak.flatten(single_e_matches_leg0, axis=-1), axis=1)
                single_electron_triggered = ak.where(trigger_fired & is_single_el, True, single_electron_triggered)
                hlt_path_fired_e[trigger.hlt_field] = ak.where(single_e_matches_leg0, trigger.id, -1)

            elif is_cross_el:
                assert trigger.n_legs == len(leg_masks) == 2
                assert abs(trigger.legs[0].pdg_id) == 11
                assert abs(trigger.legs[1].pdg_id) == 15
                
                dr_matching_e = trigger_object_matching(electrons, events.TrigObj[leg_masks[0]])
                cross_e_matches_leg0 = dr_matching_e
                
                dr_matching_tau = trigger_object_matching(taus, events.TrigObj[leg_masks[1]])
                cross_e_tau_matches_leg1 = dr_matching_tau
                
                cross_e_tau_matched = (ak.any(ak.flatten(cross_e_matches_leg0, axis=-1), axis=1) & 
                                       ak.any(ak.flatten(cross_e_tau_matches_leg1, axis=-1), axis=1))
                cross_electron_triggered = ak.where(trigger_fired & is_cross_el, True, cross_electron_triggered)
                hlt_path_fired_etau[trigger.hlt_field] = ak.where(cross_e_tau_matched, trigger.id, -1)
        
        if is_cross_tau:
            assert trigger.n_legs == len(leg_masks) == 2
            assert abs(trigger.legs[0].pdg_id) == 15
            assert abs(trigger.legs[1].pdg_id) == 15
            
            taus1 =  pair_objects.tautau.lep0
            taus2 =  pair_objects.tautau.lep1
            
            dr_matching_tau1 = trigger_object_matching(taus1, events.TrigObj[leg_masks[0]])
            dr_matching_tau2 = trigger_object_matching(taus2, events.TrigObj[leg_masks[1]])
            
            cross_tau1_matches_leg0 =  dr_matching_tau1
            cross_tau2_matches_leg1 =  dr_matching_tau2
            cross_tau_tau_matched = (ak.any(ak.flatten(cross_tau1_matches_leg0, axis=-1), axis=1) & 
                                     ak.any(ak.flatten(cross_tau2_matches_leg1, axis=-1), axis=1))
            cross_tau_triggered = ak.where(trigger_fired & is_cross_tau, True, cross_tau_triggered)
            hlt_path_fired_tau[trigger.hlt_field] = ak.where(cross_tau_tau_matched, trigger.id, -1)
    
    # Process trigger IDs for each type of lepton pair
    triggerID_e     = hlt_path_fired(events, hlt_path_fired_e)
    triggerID_etau  = hlt_path_fired(events, hlt_path_fired_etau)
    triggerID_mu    = hlt_path_fired(events, hlt_path_fired_mu)
    triggerID_mutau = hlt_path_fired(events, hlt_path_fired_mutau)
    triggerID_tau   = hlt_path_fired(events, hlt_path_fired_tau)
    
    # Generate candidate pairs based on matching triggers
    etau_channel_mask = ((triggerID_e > 0) | (triggerID_etau > 0))
    mutau_channel_mask = ((triggerID_mu > 0) | (triggerID_mutau > 0))
    tautau_channel_mask = (triggerID_tau > 0)
    
    matched_masks = DotDict.wrap({
        'etau'  : etau_channel_mask,
        'mutau' : mutau_channel_mask,
        'tautau': tautau_channel_mask
    })
    return matched_masks



