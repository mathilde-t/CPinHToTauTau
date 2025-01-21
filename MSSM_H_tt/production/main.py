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
from columnflow.production.util import attach_coffea_behavior

from MSSM_H_tt.production.weights import muon_weight, tau_weight, get_mc_weight, zpt_weight, electron_weight
from MSSM_H_tt.production.sample_split import split_dy
from MSSM_H_tt.production.generatorZ import generatorZ
from MSSM_H_tt.production.dilepton_features import hcand_fields

from MSSM_H_tt.production.aux_columns import jet_pt_def,jets_taggable,number_b_jet

np = maybe_import("numpy")
ak = maybe_import("awkward")
coffea = maybe_import("coffea")
maybe_import("coffea.nanoevents.methods.nanoaod")

# helpers
set_ak_column_f32 = functools.partial(set_ak_column, value_type=np.float32)
set_ak_column_i32 = functools.partial(set_ak_column, value_type=np.int32)

@producer(
    uses={
        attach_coffea_behavior,
        normalization_weights,
        split_dy,
        pu_weight,
        muon_weight,
        tau_weight,
        electron_weight,
        generatorZ,
        zpt_weight,
        get_mc_weight,
        hcand_fields,
        category_ids,
        jet_pt_def,
        number_b_jet,
        jets_taggable,
        "Jet.pt",
        "Jet.pt_no_jec",
        },
    produces={
        attach_coffea_behavior,
        normalization_weights,
        split_dy,
        pu_weight,
        muon_weight,
        get_mc_weight,
        tau_weight,
        electron_weight,
        generatorZ,
        zpt_weight,
        hcand_fields,
        category_ids,
        "Jet.jec_no_jec_diff",
        "Jet.number_of_jets",
        "Jet.n_jets",
        "Jet.leading_jet_pt",
        "Jet.subleading_jet_pt",
        "Jet.delta_eta_jj",
        "Jet.mjj",
        "Jet.n_jets_tag",
        
    },
)
def main(self: Producer, events: ak.Array, **kwargs) -> ak.Array:
    
    # ensure coffea behaviors are loaded
    events = self[attach_coffea_behavior](events, **kwargs)
    print("Producing Jet features...")
    events = set_ak_column_f32(events, "Jet.jec_no_jec_diff", (events.Jet.pt - events.Jet.pt_no_jec))
    events = set_ak_column_f32(events, "Jet.number_of_jets", ak.num(events.Jet))
    print("Producing Number of b-jets for categorization")
    events = self[number_b_jet](events, **kwargs)
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

        events = self[generatorZ](events, **kwargs)
        print("Z pt reweighting...")
        events = self[zpt_weight](events,**kwargs)
        print("Producing PU weights...")          
        events = self[pu_weight](events, **kwargs)
        print("Producing Muon weights...")
        events = self[muon_weight](events,do_syst = True, **kwargs)
        print("Producing Electron weights...")
        events = self[electron_weight](events,do_syst = True, **kwargs)
        print("Producing Tau weights...")
        events = self[tau_weight](events,do_syst = True, **kwargs)
        
    print("Producing jet variables for plotting...") 
    n_jets, leading_jet_pt, subleading_jet_pt, delta_eta_jj, mjj = self[jet_pt_def](events, **kwargs)
    n_jets_tag = self[jets_taggable](events, **kwargs)
    
    events = set_ak_column_f32(events, "Jet.n_jets",            n_jets           )
    events = set_ak_column_f32(events, "Jet.leading_jet_pt",    leading_jet_pt   )
    events = set_ak_column_f32(events, "Jet.subleading_jet_pt", subleading_jet_pt)
    events = set_ak_column_f32(events, "Jet.delta_eta_jj",      delta_eta_jj     )
    events = set_ak_column_f32(events, "Jet.mjj",               mjj              )
    events = set_ak_column_f32(events, "Jet.n_jets_tag",        n_jets_tag       )
    
    return events