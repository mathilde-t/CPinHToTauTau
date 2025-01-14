# coding: utf-8

"""
Column production methods related to higher-level features.
"""
import functools

from typing import Optional
from columnflow.production import Producer, producer
from columnflow.util import maybe_import
from columnflow.columnar_util import EMPTY_FLOAT, Route, set_ak_column
from columnflow.columnar_util import optional_column as optional
from httcp.util import get_lep_p4, get_ip_p4

np = maybe_import("numpy")
ak = maybe_import("awkward")
coffea = maybe_import("coffea")
maybe_import("coffea.nanoevents.methods.nanoaod")
warn = maybe_import("warnings")

# helpers
set_ak_column_f32 = functools.partial(set_ak_column, value_type=np.float32)
set_ak_column_i32 = functools.partial(set_ak_column, value_type=np.int32)

def egamma_mask(tauprod): return ((np.abs(tauprod.pdgId) == 11) + (np.abs(tauprod.pdgId) == 22))

def pion_mask(tauprod): return np.abs(tauprod.pdgId) == 211

def get_single_part(array: ak.Array, idx: int) -> ak.Array:
    return ak.firsts(array[:,idx:idx+1])

def apply_evt_mask(array: ak.Array, mask: ak.Array) -> ak.Array:
    empty_array = ak.zeros_like(array)[...][..., :0]
    return ak.where(mask, array, empty_array)
    

    
@producer(
    uses={
        "hcand_*", "tau_decay_prods_*", "event"
    },
    produces={""} 
)
def prepare_acop_vecs(self: Producer, events: ak.Array, pair_decay_ch, **kwargs):
    tau     = events.hcand_mutau.lep1 
    tauprod = events.tau_decay_prods_mutau_lep1
    muon    = events.hcand_mutau.lep0
    mask = ak.ones_like(events.hcand_mutau.lep1.pt, dtype=np.bool_)
    mask = mask & (muon.ip_sig >= 1)
    # Create masks for different channels
    if pair_decay_ch == "mu_pi":
        mask = mask & (tau.ip_sig >= 1)
        mask = mask & (tau.decayModePNet == 0) # Get only DM0 events
        mask = mask & (ak.sum(pion_mask(tauprod), axis=1) == 1) #Make sure that event conteins only one pion

    elif pair_decay_ch == "mu_rho":
        mask = mask & (tau.decayModePNet == 1) # Get only DM1 events since tau -> rho+nu -> pi^+pi^0+nu
        mask = mask & (ak.sum(pion_mask(tauprod), axis=1) == 10) # Require only one charged pion
        mask = mask & (ak.sum(egamma_mask(tauprod), axis=1) >= 1) # Require one or more electron or photon
        
    elif pair_decay_ch == "mu_a1_1pr":
        mask = mask & (tau.decayModePNet == 2) # Get only DM1 events since tau -> rho+nu -> pi^+pi^0+nu
        mask = mask & (ak.sum(pion_mask(tauprod), axis=1) == 1) # Require only one charged pion
        mask = mask & (ak.sum(egamma_mask(tauprod), axis=1) >= 1) # Require one or more electron or photon
        
    tau = apply_evt_mask(tau, ak.fill_none(ak.firsts(mask),False))
    tauprod = apply_evt_mask(tauprod, ak.fill_none(ak.firsts(mask),False))
    muon = apply_evt_mask(muon, ak.fill_none(ak.firsts(mask),False))
    # Defining a preliminary set of parameters for the function to calculate the acoplanarity angle
    # lower case variables are defined in laboratory frame
    # Get p1 and r1 that correspond to kinematic 4-vector and impact parameter vector of the muon
    p1 = get_lep_p4(muon)
    r1 = get_ip_p4(muon)
    ch1 = muon.charge
    if pair_decay_ch == "mu_pi":
        pion    = ak.mask(tauprod, pion_mask(tauprod))
        p2 = get_lep_p4(pion) # for the tau -> rho decay, p1 - is 4-vector of the charged pion and r1 is 4-vector of the neutral pion
        r2 = get_ip_p4(tau) # Create 4-vectors of tau impact parameters
        # For this channel there is no need to do the phase shift, so this arrray is filled with zeros
        do_phase_shift = ak.zeros_like(r1.energy, dtype=np.bool_)
    elif pair_decay_ch == "mu_rho":
        charged_pion = ak.mask(tauprod, pion_mask(tauprod))
        em_particles = ak.mask(tauprod, egamma_mask(tauprod))
        # for the tau -> rho decay, p1 - is 4-vector of the charged pion and r1 is 4-vector of the neutral pion
        p2 = get_lep_p4(charged_pion)
        # Create 4-vectors of tau impact parameters
        r2 = get_lep_p4(em_particles).sum()
        r2 = ak.mask(r2, r2.rho2 > 0)
        do_phase_shift = ((p2.energy - r2.energy)/(p2.energy + r2.energy)) < 0

    vecs_p4 = {}
    for var in ['p1', 'p2', 'r1', 'r2']:
        exec(f'vecs_p4["{var}"] = ak.firsts({var}, axis=1)')
    return events, vecs_p4, do_phase_shift, ch1


# def get_single_ch_tt(self: Producer, events: ak.Array,  tau_decay_channel: str,  **kwargs):
    
#     # Produce mask for different channels and tau decay modes
#     events  = ak.Array(events, behavior=coffea.nanoevents.methods.nanoaod.behavior)
#     mask    = events.channel_id == self.config_inst.get_channel('tautau').id
#     tau = []
#     tauprod = []
#     for tau_idx in range(2):
#         tau.append(get_single_part(events.hcand, tau_idx)) 
#         tauprod.append(get_single_part(events.hcandprod, tau_idx))
    
    
#     # Create masks for different channels
    
#     if tau_decay_channel == "pi_pi":
#         for tau_idx in range(2):
#             mask = mask & ak.fill_none(tau[tau_idx].decayMode == 0, False)  # Get only DM0 events
#             mask = mask & ak.fill_none(tau[tau_idx].ip_sig >= 1, False)
#             mask = mask & ak.fill_none(ak.sum(pion_mask(tauprod[tau_idx]), axis=1) == 1, False) #Make sure that event conteins only one pion

#     elif tau_decay_channel == "rho_rho":
#         for tau_idx in range(2):
#             mask = mask & ak.fill_none(tau[tau_idx].decayMode == 1, False) # Get only DM1 events since tau -> rho+nu -> pi^+pi^0+nu
#             mask = mask & ak.fill_none(ak.sum(pion_mask(tauprod[tau_idx]), axis=1) == 1, False)# Require only one charged pion
#             mask = mask & ak.fill_none(ak.sum(egamma_mask(tauprod[tau_idx]), axis=1) >= 1, False) # Require one or more electron or photon
#     hcand     = apply_evt_mask(events.hcand,     mask)
#     hcandprod = apply_evt_mask(events.hcandprod, mask)
#     # Defining a preliminary set of parameters for the function to calculate the acoplanarity angle
#     # lower case variables are defined in laboratory frame
#     # Get p1 and r1 that correspond to kinematic 4-vector and impact parameter vector of the muon
#     # Redefine taus and tauprods now from preselected events
#     tau = []
#     tauprod = []
#     for tau_idx in range(2):
#         tau.append(get_single_part(hcand, tau_idx)) 
#         tauprod.append(get_single_part(hcandprod, tau_idx))
    
#     vecs_p4 = {}
#     ch1 = tau[0].charge
#     if tau_decay_channel == "pi_pi":
#         for tau_idx in range(2):
#             charged_pion = get_single_part(ak.mask(tauprod[tau_idx], pion_mask(tauprod[tau_idx])),0)
#             vecs_p4[f"p{tau_idx+1}"] = get_lep_p4(charged_pion)
#             vecs_p4[f"r{tau_idx+1}"] = get_ip_p4(tau[tau_idx])
#         do_phase_shift = ak.zeros_like(vecs_p4["p1"].energy, dtype=np.bool_)
#     elif tau_decay_channel == "rho_rho":
#         y = []
#         for tau_idx in range(2):
#             charged_pion = get_single_part(ak.mask(tauprod[tau_idx], pion_mask(tauprod[tau_idx])),0)
#             em_particles = ak.mask(tauprod[tau_idx], egamma_mask(tauprod[tau_idx]))
#             # for the tau -> rho decay, p1 - is 4-vector of the charged pion and r1 is 4-vector of the neutral pion
#             p = get_lep_p4(charged_pion)
#             r = get_lep_p4(em_particles).sum()
#             r = ak.mask(r, r.rho2 > 0)
#             vecs_p4[f"r{tau_idx+1}"] = p
#             vecs_p4[f"p{tau_idx+1}"] = r
#             y.append((p.energy - r.energy)/(p.energy + r.energy))
#         do_phase_shift = (y[0] * y[1]) < 0

#     return events, vecs_p4, do_phase_shift, ch1


def make_boost(vecs_p4, boostvec=None):
    # Check if  boost vector is provided, othervise calculated it by getting p1 and p2 from the vecs_p4 dict
    if boostvec == None:
        boostvec_ = vecs_p4['p1'] + vecs_p4['p2']
    else:
        boostvec_ = boostvec
    # Create a dictionary to store boosted variables (they are defined with upper case names)
    zmf_vars = {}
    for var in vecs_p4.keys():
        exec(f'zmf_vars["{var.upper()}"] = vecs_p4["{var}"].boost(boostvec_.boostvec.negative())')
    return zmf_vars


def get_acop_angle(vecs_p4, do_phase_shift):
    # Create 4 3-vectors from the vecs_p4 dict
    P1 = vecs_p4['P1'].pvec.unit
    P2 = vecs_p4['P2'].pvec.unit
    R1 = vecs_p4['R1'].pvec.unit
    R2 = vecs_p4['R2'].pvec.unit
    R1_tan = R1.add((P1.multiply(R1.dot(P1))).negative())
    R2_tan = R2.add((P2.multiply(R2.dot(P2))).negative())

    O = P2.dot(R1_tan.cross(R2_tan))

    raw_phi = np.acos((R1_tan.unit).dot(R2_tan.unit))
    phi_cp = ak.where(O > 0, raw_phi, 2 * np.pi - raw_phi)
    phi_cp = ak.where(do_phase_shift, phi_cp + np.pi, phi_cp)
    # Map  angles into [0,2pi) interval
    phi_cp = ak.where(phi_cp > 2.*np.pi, phi_cp - 2. * np.pi, phi_cp)
    phi_cp = ak.where(phi_cp < 0, phi_cp + 2. * np.pi, phi_cp)
    return phi_cp


def produce_dy_alpha(
        vecs_p4,
        ch1) -> ak.Array:
    empty_p2 = ak.zeros_like(vecs_p4['p2'].pvec.unit)[...][..., :0]
    empty_r2 = ak.zeros_like(vecs_p4['r2'].pvec.unit)[...][..., :0]
    mask = ak.fill_none(ak.firsts(ch1 < 0), False)
    P_ = ak.where(mask, vecs_p4['p1'].pvec.unit, vecs_p4['p2'].pvec.unit)
    R_ = ak.where(mask, vecs_p4['r1'].pvec.unit, vecs_p4['r2'].pvec.unit)
    z = ak.zeros_like(P_)
    z['z'] = 1
    vec1 = (z.cross(P_)).unit
    vec2 = (R_.cross(P_)).unit
    dy_alpha = np.arccos(np.absolute(vec1.dot(vec2)))
    dy_alpha = ak.where(np.isnan(dy_alpha),
                        ak.ones_like(dy_alpha)*EMPTY_FLOAT,
                        dy_alpha)
    return dy_alpha


channels = ['mu_pi','mu_rho']


@producer(
    uses={
        prepare_acop_vecs,
        "hcand_*",

        optional("GenTau.*"), optional("GenTauProd.*")},
    produces={
        f"phi_cp_{the_ch}" for the_ch in channels
    } | {
        optional(f"phi_cp_{the_ch}_reg1") for the_ch in channels
    } | {
        optional(f"phi_cp_{the_ch}_reg2") for the_ch in channels
    } | {
        f"phi_cp_{the_ch}_2bin" for the_ch in channels
    } | {
        optional(f"phi_cp_{the_ch}_reg1_2bin") for the_ch in channels
    } | {
        optional(f"phi_cp_{the_ch}_reg2_2bin") for the_ch in channels
    } | {
        optional(f"alpha_{the_ch}") for the_ch in channels
    } | {prepare_acop_vecs}

)
def phi_cp(
        self: Producer,
        events: ak.Array,
        **kwargs
) -> ak.Array:

    alpha = np.ones_like(events.event)*EMPTY_FLOAT
    for the_ch in channels:
        print(f"Calculating phi_cp for {the_ch}")
        if 'mu' in the_ch: events, vecs_p4, do_phase_shift, ch1 = self[prepare_acop_vecs](events, the_ch, **kwargs)
        zmf_vecs_p4 = make_boost(vecs_p4)
        phi_cp = get_acop_angle(zmf_vecs_p4, do_phase_shift)
        # alpha_per_ch = produce_dy_alpha(vecs_p4, ch1)
        # alpha = ak.where(ak.num(alpha_per_ch) > 0,
        #                  ak.fill_none(ak.firsts(alpha_per_ch), EMPTY_FLOAT),
        #                  alpha)
        # reg1_mask = ak.fill_none(ak.firsts(alpha_per_ch >= np.pi/4.), False)
        # reg2_mask = ak.fill_none(
        #     ak.firsts((alpha_per_ch < np.pi/4.) & (alpha_per_ch >= 0.)), False)

        phi_cp_2bin = (phi_cp + np.pi/2.) % (2*np.pi)
        events = set_ak_column_f32(
            events, f"phi_cp_{the_ch}_2bin",  ak.fill_none(ak.firsts(phi_cp_2bin),EMPTY_FLOAT))
        # for the_reg in ['reg1', 'reg2']:
        #     var = ak.where(eval(f'{the_reg}_mask'), phi_cp, EMPTY_FLOAT)
        #     var_2bin = ak.where(
        #         eval(f'{the_reg}_mask'), phi_cp_2bin, EMPTY_FLOAT)

        #     events = set_ak_column_f32(
        #         events, f"phi_cp_{the_ch}_{the_reg}",  var)
        #     events = set_ak_column_f32(events, f"phi_cp_{the_ch}_{the_reg}_2bin",  var_2bin)

        # events = set_ak_column_f32(
        #     events, f'alpha_{the_ch}', ak.fill_none(alpha, EMPTY_FLOAT))
        events = set_ak_column_f32(
            events, f"phi_cp_{the_ch}",  ak.fill_none(ak.firsts(phi_cp),EMPTY_FLOAT))
    return events
