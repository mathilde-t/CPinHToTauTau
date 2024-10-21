"""
Produce channel_id column. This function is called in the main selector
"""

from columnflow.production import Producer, producer
from columnflow.selection import Selector, SelectionResult, selector
from columnflow.columnar_util import set_ak_column
from columnflow.util import maybe_import, DotDict


np = maybe_import("numpy")
ak = maybe_import("awkward")

import functools
set_ak_column_f32 = functools.partial(set_ak_column, value_type=np.float32)
set_ak_column_i32 = functools.partial(set_ak_column, value_type=np.int32)


@producer(
    produces={
        "channel_id",
    },
    exposed=False,
)
def channel_id(
        self: Producer,
        events: ak.Array,
        **kwargs
) -> ak.Array:
    #Create a template array filled with zeros
    channel_id  = ak.zeros_like(ak.local_index(events.event), dtype = np.uint8)
    
    #Create a mask array to check for channel orthogonality 
    and_mask = ak.ones_like(ak.local_index(events.event), dtype = np.bool_)
    
    for channel in self.config_inst.channels.names():
        the_mask = ak.num(events[f'hcand_{channel}'].lep0, axis=1) > 0
        and_mask = and_mask & the_mask
        channel_id = ak.where(the_mask,
                              self.config_inst.get_channel(channel).id,
                              channel_id)

    if ak.any(and_mask):
        raise TypeError('Found events belonging to multiple channels!')
    else:
        channel_id = ak.values_astype(channel_id, np.uint8)
        events = set_ak_column(events, "channel_id", channel_id)

    return events





@producer(
    uses={ f"Jet.{var}" for var in 
        [
            "pt", "eta", "phi", "mass",
            "jetId", "btagDeepFlavB"
        ]} | {"hcand_*"},
    produces={"is_b_vetoed"},
    exposed=False,
)
def jet_veto(
        self: Producer,
        events: ak.Array,
        **kwargs
) -> ak.Array:
    """
    This function that produces 'is_b_vetoed' column to be used in categorisation
    """
    year = self.config_inst.x.year
    btag_wp = self.config_inst.x.btag_working_points[year].deepjet.medium
    # nominal jet selection
    jet_selections = {
        "jet_pt_20"               : events.Jet.pt > 20.0,
        "jet_eta_2.5"             : abs(events.Jet.eta) < 2.5,
        "jet_id"                  : events.Jet.jetId & 0b10 == 0b10, 
        "btag_wp_medium"          : events.Jet.btagDeepFlavB >= btag_wp
    }
    jet_obj_mask = ak.ones_like(ak.local_index(events.Jet.pt, axis=1), dtype=np.bool_)
    for the_sel in jet_selections.values():
        jet_obj_mask = jet_obj_mask & the_sel
    presel_jet = ak.drop_none(ak.mask(events.Jet, jet_obj_mask))
    
    leg_masks = []
    event_mask = ak.zeros_like(events.event, dtype=np.bool_)
    for the_ch in self.config_inst.channels.names(): 
        hcand = events[f'hcand_{the_ch}']
        for idx in range(2): #iteate over taus
            lep = ak.firsts(hcand[f'lep{idx}'])
            jet_tau_pairs = ak.cartesian([presel_jet,lep], axis=1)
            jet_br, lep_br = ak.unzip(jet_tau_pairs)
            delta_r = jet_br.delta_r(lep_br)
            event_mask = event_mask | ak.any(delta_r < 0.4, axis=1)  #make joint mask for first and second tau.
    events = set_ak_column(events, "is_b_vetoed", event_mask)
    return events




def is_pion(prods): return ((np.abs(prods.pdgId) == 211))

def is_photon(prods): return prods.pdgId == 22

def has_one_pion(prods): return (ak.sum(is_pion(prods),   axis=1) == 1)

def has_three_pions(prods): return (ak.sum(is_pion(prods),   axis=1) == 3)

def has_photons(prods): return (ak.sum(is_photon(prods), axis=1) > 0)

def has_no_photons(prods): return (ak.sum(is_photon(prods), axis=1) == 0)


@producer(
    uses={
        "TauProd.pdgId",
    },
    produces={
        "TauProd.mass", "TauProd.charge",
    },
    exposed=False,
)
def assign_tauprod_mass_charge(
        self: Selector,
        events: ak.Array,
        **kwargs
) -> ak.Array:
    
    #https://pdg.lbl.gov/2023/listings/particle_properties.html
    
    part_dict=DotDict.wrap({
        "pion_pm"   : {'mass'    : 0.13957039, #GeV
                       'pdg_id'  : 211},
        "pion_0"    : {'mass'    : 0.1349768, #GeV
                       'pdg_id'  : 111},
        "gamma"     : {'mass'    : 0.0, #GeV
                       'pdg_id'  : 22},
        "kaon_pm"   : {'mass'    : 0.493677, #GeV
                       'pdg_id'  : 321},
        "ele_pm"    : {'mass'    : 0.00051099895, #GeV
                       'pdg_id'  : 11},
        "muon_pm"   : {'mass'    : 0.1056583755, #GeV
                       'pdg_id'  : 13},
        })
    mass = ak.zeros_like(ak.local_index(events.TauProd.pdgId), dtype=np.float32)
    charge = ak.zeros_like(ak.local_index(events.TauProd.pdgId), dtype=np.int32)
    for part in part_dict:
       
        prod_id = events.TauProd.pdgId
        mass = ak.where(np.abs(prod_id)==part_dict[part].pdg_id,
                        part_dict[part].mass,
                        mass)
        if '_pm' in part:
            charge = ak.where(np.abs(prod_id)==part_dict[part].pdg_id,
                              np.sign(prod_id),
                              charge)
    
    events = set_ak_column_f32(events, "TauProd.mass", mass)
    events = set_ak_column_i32(events, "TauProd.charge", charge)
    return events



@producer(
    uses={
        "TauProd.*", assign_tauprod_mass_charge,
    },
    produces={
        "tau_decay_prods_*",
        assign_tauprod_mass_charge,
    },
    exposed=False,
)
def add_tau_prods(
        self: Producer,
        events: ak.Array,
        **kwargs
) -> tuple[ak.Array, SelectionResult]:
    
    events   = self[assign_tauprod_mass_charge](events)
    channels = self.config_inst.channels.names()
    ch_objects = self.config_inst.x.ch_objects
    mask = ak.zeros_like(events.event, dtype=np.bool_)
    for ch_str in channels: 
        hcand = events[f'hcand_{ch_str}']
        tau_decay_prods_dict = {}
        for lep_str in hcand.fields:
            if ch_objects[ch_str][lep_str] == 'Tau':
                tau = hcand[lep_str]
                tauprods = events.TauProd
                tau_idx = tau.rawIdx
                tauprods_mask = tauprods.tauIdx
                idx_pairs = ak.cartesian([tau_idx,tauprods.tauIdx], axis=1, nested=True)
                tau_idx_b, tau_prod_idx_b = ak.unzip(idx_pairs)
                tau2prod_match_mask = ak.fill_none(tau_idx_b == tau_prod_idx_b, False)
                matched_tau_prods = tauprods[ak.flatten(tau2prod_match_mask,axis=2)]
                #DM0
                tau = ak.firsts(tau)
                mask = mask | ak.fill_none(tau.decayMode==0, False) & has_one_pion(matched_tau_prods)
                #DM1
                mask = mask | ak.fill_none(tau.decayMode==1, False) & has_one_pion(matched_tau_prods) & has_photons(matched_tau_prods)
                #DM10
                mask = mask | ak.fill_none(tau.decayMode==10, False) & has_three_pions(matched_tau_prods) 
                #DM11
                mask = mask | ak.fill_none(tau.decayMode==11, False) & has_three_pions(matched_tau_prods) & has_photons(matched_tau_prods)
                events = set_ak_column(events, f'tau_decay_prods_{ch_str}_{lep_str}',  matched_tau_prods)
            else:
                pass
        
    return events, SelectionResult(
        steps={
            "decay_prods_are_ok": mask,
        },
    )

