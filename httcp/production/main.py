# coding: utf-8

"""
Column production methods related to higher-level features.
"""
import functools

from typing import Optional
from columnflow.production import Producer, producer
from columnflow.production.categories import category_ids
from columnflow.production.normalization import normalization_weights
from columnflow.production.cms.pileup import pu_weight
from columnflow.production.cms.seeds import deterministic_seeds
from columnflow.selection.util import create_collections_from_masks
from columnflow.util import maybe_import
from columnflow.columnar_util import EMPTY_FLOAT, Route, set_ak_column
from columnflow.columnar_util import optional_column as optional

#from httcp.production.PhiCPNeutralPion import PhiCPNPMethod
from httcp.production.ReArrangeHcandProds import reArrangeDecayProducts, reArrangeGenDecayProducts
from httcp.production.PhiCP_Producer import ProduceDetPhiCP, ProduceGenPhiCP

from httcp.production.weights import muon_weight, tau_weight, get_mc_weight,tauspinner_weight
from httcp.production.sample_split import split_dy
from httcp.production.dilepton_features import hcand_mass, mT, rel_charge, hcand_mt
from httcp.production.phi_cp import phi_cp
from httcp.util import get_lep_p4

np = maybe_import("numpy")
ak = maybe_import("awkward")
coffea = maybe_import("coffea")
maybe_import("coffea.nanoevents.methods.nanoaod")

# helpers
set_ak_column_f32 = functools.partial(set_ak_column, value_type=np.float32)
set_ak_column_i32 = functools.partial(set_ak_column, value_type=np.int32)


@producer(
    uses={
        'hcand_*', 'PuppiMET*'
    },
    produces={
        'hcand_*'
    },
)
def hcand_fields(
        self: Producer, 
        events: ak.Array,
        **kwargs
) -> ak.Array:
    channels = self.config_inst.channels.names()
    ch_objects = self.config_inst.x.ch_objects
    for ch_str in channels:
        hcand = events[f'hcand_{ch_str}']
        p4 = {}
        for the_lep in hcand.fields: p4[the_lep] = get_lep_p4(hcand[the_lep]) 
        
        mass = (p4['lep0'] + p4['lep1']).mass
        hcand['mass'] = ak.where(mass > 0, mass , EMPTY_FLOAT)
        
        delta_r = ak.flatten(p4['lep0'].metric_table(hcand.lep1), axis=2)
        hcand['delta_r'] = ak.where(delta_r > 0, delta_r , EMPTY_FLOAT)
        hcand['rel_charge'] = hcand.lep0.charge * hcand.lep1.charge
        if ch_str !=' tautau':
            mt = hcand_mt(p4['lep0'], events.PuppiMET)
            hcand['mt'] = ak.where(mt > 0, mt , EMPTY_FLOAT)   
        events = set_ak_column(events, f'hcand_{ch_str}', hcand) 
    # events, P4_dict     = self[reArrangeDecayProducts](events)
    # events              = self[ProduceDetPhiCP](events, P4_dict)

    # if "is_signal" in list(self.dataset_inst.aux.keys()):
    #     if self.dataset_inst.aux["is_signal"]:
    #         events, P4_gen_dict = self[reArrangeGenDecayProducts](events)
    #         events = self[ProduceGenPhiCP](events, P4_gen_dict)
    return events


@producer(
    uses={
        normalization_weights,
        split_dy,
        pu_weight,
        muon_weight,
        tau_weight,
        get_mc_weight,
        hcand_fields,
        tauspinner_weight,
        phi_cp,
        category_ids,
    },
    produces={
        normalization_weights,
        split_dy,
        pu_weight,
        muon_weight,
        get_mc_weight,
        tau_weight,
        hcand_fields,
        tauspinner_weight,
        phi_cp,
        category_ids,
    },
)
def main(self: Producer, events: ak.Array, **kwargs) -> ak.Array:
    print("Producing Hcand features...")
    events = self[hcand_fields](events, **kwargs) 
    events = self[category_ids](events, **kwargs)
    if self.dataset_inst.is_mc:
        events = self[get_mc_weight](events, **kwargs)
        print("Producing Normalization weights...")
        events = self[normalization_weights](events, **kwargs)
        processes = self.dataset_inst.processes.names()
        if ak.any(['dy' in proc for proc in processes]):
            print("Splitting Drell-Yan dataset...")
            events = self[split_dy](events,**kwargs)
        print("Producing PU weights...")  
        events = self[pu_weight](events, **kwargs)
        print("Producing Muon weights...")
        events = self[muon_weight](events,do_syst = True, **kwargs)
        print("Producing Tau weights...")
        events = self[tau_weight](events,do_syst = True, **kwargs)
        print("Producing Tauspinner weights...")
        events = self[tauspinner_weight](events, **kwargs)
    print("Producing phi_cp...") 
    events = self[phi_cp](events, **kwargs) 
    return events