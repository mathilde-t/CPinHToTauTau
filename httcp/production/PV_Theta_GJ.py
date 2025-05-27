"""
Adaption of the Polarimetic Vector Estimation via use of the Goldfried-Jackson Angle
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
            },
    produces={
        # new columns
        # "theta_GJ_4vector",  # pour l'angle 4 vectoriel
        "theta_GJ_scalar",  # pour l'angle basé sur énergie

        "E_tau" ,"E_vis" ,
        "m_tau" ,"m_vis" ,
        "p_tau" ,"p_vis" ,
    },
)

def compute_theta_GJ(
        self: Producer, 
        events: ak.Array,
        **kwargs
) -> ak.Array:
    """
    Calculates the Gottfried-Jackson angle θ_GJ between the a1 meson (in the tau rest frame)
    and the tau momentum (in the Higgs rest frame). This angle carries spin information.

    Args:
        events: awkward Array with NanoAOD behavior and hcand/tau info

    Returns:
        events: awkward Array with new column "theta_GJ" [float32]
        using : cos theta_GJ = ( 2E_tau E_vis - m_tau^2 - m_vis^2 ) / (2 p_tau p_vis )
    """


    # Convert events to awkward array with nanoaod behavior
    events = ak.Array(events, behavior=coffea.nanoevents.methods.nanoaod.behavior)


    hcand_ = ak.with_name(events.hcand, "PtEtaPhiMLorentzVector")
    hcand1 = hcand_[:, 0:1] # leading tau and the muon in the channel : tautau -> mu a1
    hcand2 = hcand_[:, 1:2]

    E_tau = hcand1.energy
    E_vis = hcand2.energy
    m_tau = hcand1.mass
    m_vis = hcand2.mass
    p_tau = hcand1.pt
    p_vis = hcand2.pt

    # Calculate cos(theta_GJ) using the scalar formula
    cos_theta_GJ_scalar = (2 * E_tau * E_vis - m_tau**2 - m_vis**2) / (2 * p_tau * p_vis)
    cos_theta_GJ_scalar = ak.clip(cos_theta_GJ_scalar, -1.0, 1.0)
    theta_GJ_scalar = ak.arccos(cos_theta_GJ_scalar)
    logger.debug("theta_GJ_scalar (mean): %s", ak.mean(theta_GJ_scalar))

    # Store scalar result
    set_ak_column_f32(events, "theta_GJ_scalar", theta_GJ_scalar)
    set_ak_column_f32(events, "E_tau", E_tau)
    set_ak_column_f32(events, "E_vis", E_vis)
    set_ak_column_f32(events, "m_tau", m_tau)
    set_ak_column_f32(events, "m_vis", m_vis)
    set_ak_column_f32(events, "p_tau", p_tau)
    set_ak_column_f32(events, "p_vis", p_vis)


    return events



    # def get_beta(vec):
#     """Return 3-velocity β from LorentzVector"""
#     return ak.zip(
#         {
#             "x": vec.x / vec.t,
#             "y": vec.y / vec.t,
#             "z": vec.z / vec.t,
#         },
#         with_name="ThreeVector",
#         behavior=coffea.nanoevents.methods.vector.behavior,
#     )

# def dot_product(v1, v2):
#     return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z

# def norm(v):
#     return np.sqrt(v.x**2 + v.y**2 + v.z**2)

# def unit_vector(v):
#     n = norm(v)
#     return ak.zip(
#         {
#             "x": v.x / n,
#             "y": v.y / n,
#             "z": v.z / n,
#         },
#         with_name="ThreeVector",
#         behavior=coffea.nanoevents.methods.vector.behavior,
#     )


    # # Use rearranged decay product info
    # events, decay_products = self.reArrangeDecayProducts(events)

    # # Get tau 4-vectors
    # tau1 = decay_products["p4h1"]
    # tau2 = decay_products["p4h2"]
    # higgs_cand = tau1 + tau2
    # hmass = higgs_cand.mass

    # # Build the a1 vector — assuming it decays into 3 charged pions from tau2
    # #a1 = ak.sum(decay_products["p4h2pi"], axis=1)
    # valid = ak.num(decay_products["p4h2pi"]) >= 3
    # a1 = ak.sum(decay_products["p4h2pi"], axis=1)
    # a1 = ak.where(valid, a1, EMPTY_FLOAT)

    # # Boost a1 into tau2 rest frame
    # beta_tau2 = get_beta(tau2)
    # a1_in_tau2_RF = a1.boost(-beta_tau2)

    # # Boost tau2 into Higgs rest frame
    # beta_higgs = get_beta(higgs_cand)
    # tau2_in_H_RF = tau2.boost(-beta_higgs)

    # # Extract direction unit vectors
    # a1_dir = unit_vector(a1_in_tau2_RF)
    # tau2_dir = unit_vector(tau2_in_H_RF)

    # # Compute angle between vectors
    # cos_theta_GJ = dot_product(a1_dir, tau2_dir)
    # cos_theta_GJ = ak.clip(cos_theta_GJ, -1.0, 1.0)
    # theta_GJ = ak.arccos(cos_theta_GJ)
    # logger.debug("theta_GJ (mean): %s", ak.mean(theta_GJ))

    # # Store result
    # set_ak_column_f32(events, "theta_GJ", theta_GJ)