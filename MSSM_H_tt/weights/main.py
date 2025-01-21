# coding: utf-8

"""
Example event weight producer.
"""

from columnflow.weight import WeightProducer, weight_producer
from columnflow.util import maybe_import
from columnflow.config_util import get_shifts_from_sources
from columnflow.columnar_util import Route

ak = maybe_import("awkward")
np = maybe_import("numpy")


@weight_producer(
    # both used columns and dependent shifts are defined in init below
    # only run on mc
    mc_only=True,
)
def main(self: WeightProducer, events: ak.Array, **kwargs) -> ak.Array:
    # build the full event weight
    weight = ak.Array(np.ones(len(events), dtype=np.float32))
    for column in self.weight_columns:
        weight = weight * Route(column).apply(events)
    processes = self.dataset_inst.processes.names()
    process_id = events.process_id
    Z_ee_weight = 1
    if ak.any(['dy' in proc for proc in processes]) and self.dataset_inst.campaign.x.tag == 'postEE':
        print("Applying an ad hoc weight to Z->ee...")
        weight = ak.where(process_id==51001,weight*Z_ee_weight,weight)
        print("The ad hoc weight to Z->ee was applied succefully...")
        
    return events, weight


@main.init
def main_init(self: WeightProducer) -> None:
    # store column names referring to weights to multiply
    self.weight_columns = {
        "normalization_weight",
        "mc_weight",
        "pu_weight",
        "muon_weight_nom",
        "tau_weight_nom",
        "electron_weight_nom",
        "zpt_weight"
    }
    self.uses |= self.weight_columns
    
    # declare shifts that the produced event weight depends on
    shift_sources = {}
    self.shifts |= set(get_shifts_from_sources(self.config_inst, *shift_sources))