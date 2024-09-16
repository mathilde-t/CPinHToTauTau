"""
Produce channel_id column. This function is called in the main selector
"""

from columnflow.production import Producer, producer
from columnflow.columnar_util import set_ak_column
from columnflow.util import maybe_import

np = maybe_import("numpy")
ak = maybe_import("awkward")


# @selector(
#     #uses={#
#     #
#     #   },
#     produces={
#         "channel_id",
#     },
#     exposed=False,
# )
# def get_categories(
#         self: Selector,
#         events: ak.Array,
#         trigger_results: SelectionResult,
#         etau_pair_indices: ak.Array,
#         mutau_pair_indices: ak.Array,
#         tautau_pair_indices: ak.Array,
#         **kwargs
# ) -> tuple[ak.Array, SelectionResult]:
#     # get channels from the config
#     ch_etau   = self.config_inst.get_channel("etau")
#     ch_mutau  = self.config_inst.get_channel("mutau")
#     ch_tautau = self.config_inst.get_channel("tautau")

#     false_mask       = (abs(events.event) < 0)
#     single_triggered = false_mask
#     cross_triggered  = false_mask
#     empty_indices    = ak.zeros_like(1 * events.event, dtype=np.uint16)[..., None][..., :0]

#     channel_selections = {
#         "cat_is_etau"           : [ch_etau.id, 
#                                    ((ak.num(etau_pair_indices, axis=1) == 2) 
#                                     & (ak.num(mutau_pair_indices, axis=1) == 0) 
#                                     & (ak.num(tautau_pair_indices, axis=1) == 0))],
#         "cat_is_mutau"          : [ch_mutau.id, 
#                                    ((ak.num(etau_pair_indices, axis=1) == 0) 
#                                     & (ak.num(mutau_pair_indices, axis=1) == 2) 
#                                     & (ak.num(tautau_pair_indices, axis=1) == 0))],
#         "cat_is_tautau"         : [ch_tautau.id, 
#                                    ((ak.num(etau_pair_indices, axis=1) == 0) 
#                                     & (ak.num(mutau_pair_indices, axis=1) == 0) 
#                                     & (ak.num(tautau_pair_indices, axis=1) == 2))],
#         "cat_is_etau_mutau"     : [ch_etau.id + ch_mutau.id, 
#                                    ((ak.num(etau_pair_indices, axis=1) == 2) 
#                                     & (ak.num(mutau_pair_indices, axis=1) == 2) 
#                                     & (ak.num(tautau_pair_indices, axis=1) == 0))],
#         "cat_is_etau_tautau"    : [ch_etau.id + ch_tautau.id, 
#                                    ((ak.num(etau_pair_indices, axis=1) == 2) 
#                                     & (ak.num(mutau_pair_indices, axis=1) == 0) 
#                                     & (ak.num(tautau_pair_indices, axis=1) == 2))], 
#         "cat_is_mutau_tautau"   : [ch_mutau.id + ch_tautau.id, 
#                                    ((ak.num(etau_pair_indices, axis=1) == 0) 
#                                     & (ak.num(mutau_pair_indices, axis=1) == 2) 
#                                     & (ak.num(tautau_pair_indices, axis=1) == 2))], 
#     }

#     selection_steps  = {}
#     channel_id       = np.uint8(1) * false_mask
#     for key, val in channel_selections.items():
#         selection_steps[key] = val[1]
#         channel_id = ak.where(val[1], val[0], channel_id)
        
#     # some final type conversions
#     channel_id = ak.values_astype(channel_id, np.uint8)
#     events = set_ak_column(events, "channel_id", channel_id)

#     return events, SelectionResult(aux=selection_steps)



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
