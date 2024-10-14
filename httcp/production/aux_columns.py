"""
Produce channel_id column. This function is called in the main selector
"""

from columnflow.production import Producer, producer
from columnflow.columnar_util import set_ak_column
from columnflow.util import maybe_import

np = maybe_import("numpy")
ak = maybe_import("awkward")

@producer(
    produces={
        "channel_id",
    },
    exposed=False,
)
def channel_id(
        self: Producer,
        events: ak.Array,
        etau_mask: ak.Array,
        mutau_mask: ak.Array,
        tautau_mask: ak.Array,
        **kwargs
) -> ak.Array:
    #Create a template array filled with zeros
    channel_id  = ak.zeros_like(ak.local_index(events.event), dtype = np.uint8)
    
    #Create a mask array to check for channel orthogonality 
    and_mask = ak.ones_like(ak.local_index(events.event), dtype = np.bool_)
    for channel in ['etau', 'mutau','tautau']:
        the_mask = eval(f'{channel}_mask')
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
        ]} | {"hcand.*"},
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
    for idx in range(2): #iteate over taus
        tau = ak.firsts(events.hcand[:,idx:idx+1])
        jet_tau_pairs = ak.cartesian([presel_jet,tau], axis=1)
        jet_br, tau_br = ak.unzip(jet_tau_pairs)
        delta_r = jet_br.delta_r(tau_br)
        leg_masks.append(ak.any(delta_r < 0.4, axis=1))
    event_mask = leg_masks[0] | leg_masks[1] #make joint mask for first and second tau.
    events = set_ak_column(events, "is_b_vetoed", event_mask)
    return events 