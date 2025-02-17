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

np = maybe_import("numpy")
ak = maybe_import("awkward")
coffea = maybe_import("coffea")
maybe_import("coffea.nanoevents.methods.nanoaod")

@producer(
        uses={
            # nano columns | input for compute_fastmtt function
            "hcand.*", "channel_id",                    #includes "hcand.decayMode"
            "PuppiMET.pt", "PuppiMET.phi", "PuppiMET.covXX", "PuppiMET.covXY", "PuppiMET.covYY",

        },
        produces={
            # new columns | output of the compute_fastmtt function to use in the cf framework
            "hcand_invm_fastMTT_IC", "hcand.pt_fastMTT_IC", 
            "hcand.pt1_fastMTT_IC", "hcand.pt2_fastMTT_IC"
        }
)

def apply_fastMTT_IC(
        self: Producer,
        events: ak.Array,
        **kwargs
) -> ak.Array:

    """
    fastMTT implementation for the TIDAL project from the Imperial College London HEP group.
    https://github.com/Ksavva1021/TIDAL/blob/main/Tools/FastMTT/fastmtt.py
    Adapted to the IPHCtau coloumflow framework by the IPHCtau HiggsCP to tautau group.
    """

    events = ak.Array(events, behavior=coffea.nanoevents.methods.nanoaod.behavior)
    hcand_ = ak.with_name(events.hcand, "PtEtaPhiMLorentzVector")
    hmass  = (hcand_[:,:1] + hcand_[:,1:2]).mass
    
    etau_id   = self.config_inst.get_channel("etau").id
    mutau_id  = self.config_inst.get_channel("mutau").id
    tautau_id = self.config_inst.get_channel("tautau").id
    
    decay_type_0 = ak.zeros_like(hmass) #etau
    decay_type_1 = ak.ones_like(hmass) #mutau
    decay_type_2 = 2 * decay_type_1 #tautau

    N = len(events)
    pt_1    = ak.to_numpy(events.hcand.pt[:,0:1])
    eta_1   = ak.to_numpy(events.hcand.eta[:,0:1])
    phi_1   = ak.to_numpy(events.hcand.phi[:,0:1])
    mass1  = ak.to_numpy(events.hcand.mass[:,0:1])
    dm1    = events.hcand.decayMode[:,0:1]
    dm1_dummy = ak.values_astype(-1 * ak.ones_like(dm1), np.int32)
    dm1    = ak.to_numpy(ak.where(events.channel_id < 4, dm1_dummy, dm1))
    type1  = ak.to_numpy(ak.where(events.channel_id == etau_id,
                                  decay_type_0,
                                  ak.where(events.channel_id == mutau_id,
                                           decay_type_1,
                                           decay_type_2)
                                  )
                         )
    
    pt_2    = ak.to_numpy(events.hcand.pt[:,1:2])
    eta_2   = ak.to_numpy(events.hcand.eta[:,1:2])
    phi_2   = ak.to_numpy(events.hcand.phi[:,1:2])
    mass2  = ak.to_numpy(events.hcand.mass[:,1:2])
    dm2    = ak.to_numpy(events.hcand.decayMode[:,1:2])
    type2  = ak.to_numpy(ak.where(events.channel_id == tautau_id, decay_type_2, decay_type_0)) #0 for et and mt, 2 for tt IC : https://github.com/Ksavva1021/TIDAL/blob/main/Tools/FastMTT/fastmtt.py#L40-L48
    
    metpt  = ak.to_numpy(events.PuppiMET.pt[:,None])
    metphi = ak.to_numpy(events.PuppiMET.phi[:,None])
    metcov_xx = ak.to_numpy(events.PuppiMET.covXX[:,None])
    metcov_xy = ak.to_numpy(events.PuppiMET.covXY[:,None])
    metcov_yx = metcov_xy
    metcov_yy = ak.to_numpy(events.PuppiMET.covYY[:,None])

    met_x = metpt * np.cos(metphi)
    met_y = metpt * np.sin(metphi)

    # initialize global parameters
    m_ele = 0.51099895/1000  #MeV -> GeV 
    m_muon = 105.6583755/1000 #MeV -> GeV 
    m_tau = 1776.85/1000 #MeV -> GeV 
    m_pion = 139.5/1000 #MeV -> GeV
    delta = 1/1.15
    reg_order = 6
    constrain = True
    constrain_setting = "Window"
    constrain_window = np.array([123, 127])

    fastmttMass_values, fastmttPt_values, fastmttPt1_values, fastmttPt2_values = compute_fast_MTT(
        N,
        pt_1,eta_1,phi_1,mass1,
        pt_2,eta_2,phi_2,mass2,
        met_x,met_y,metcov_xx,metcov_xy,metcov_yx,metcov_yy,
        decay_type_1,decay_type_2,
        m_ele,
        m_muon,
        m_tau,
        m_pion,
        delta,
        reg_order,
        constrain,
        constrain_setting,
        constrain_window
    )

    events = set_ak_column(events, "hcand_invm_fastMTT_IC", fastmttMass_values)
    events = set_ak_column(events, "hcand.pt_fastMTT_IC", fastmttPt_values)
    events = set_ak_column(events, "hcand.pt1_fastMTT_IC", fastmttPt1_values)
    events = set_ak_column(events, "hcand.pt2_fastMTT_IC", fastmttPt2_values)

    #from IPython import embed; embed()

    return events

def compute_fast_MTT(
        N,
        pt_1,eta_1,phi_1,mass1,
        pt_2,eta_2,phi_2,mass2,
        met_x,met_y,metcov_xx,metcov_xy,metcov_yx,metcov_yy,
        decay_type_1,decay_type_2,
        m_ele,
        m_muon,
        m_tau,
        m_pion,
        delta,
        reg_order,
        constrain,
        constrain_setting,
        constrain_window):
    fastmttMass_values = np.zeros(N, dtype=np.float32)
    fastmttPt_values = np.zeros(N, dtype=np.float32)
    fastmttPt1_values = np.zeros(N, dtype=np.float32)
    fastmttPt2_values = np.zeros(N, dtype=np.float32)

    mass_dict = {0: m_ele, 1: m_muon, 2: m_tau}


    for i in range(70):#N):

        print("i : ", i)

        #from IPython import embed; embed()

        # grab the correct masses based on tau decay type
        # tau decay_type: 0 ==> leptonic to electron,
        #                 1 ==> leptonic to muon,
        #                 2 ==> leptonic to hadronic

        decay_type_1_i = int(ak.firsts(decay_type_1)[i]) # decay_type_1 has a nested structure : <Array [[1], [1], [1], [1], ..., [1], [1], [1]] type='5000 * var * float64'>
        decay_type_2_i = int(ak.firsts(decay_type_2)[i])

        if (decay_type_1_i != 2):
            m1 = mass_dict[decay_type_1_i]
        else:
            m1 = float(mass1[i])
        if (decay_type_2_i != 2):
            m2 = mass_dict[decay_type_2_i]
        else:
            m2 = float(mass2[i])

        # store visible masses
        m_vis_1 = m1
        m_vis_2 = m2

        # determine minimum and maximum possible masses
        m_vis_min_1, m_vis_max_1 = 0, 0
        m_vis_min_2, m_vis_max_2 = 0, 0
        if (decay_type_1_i == 0):
            m_vis_min_1, m_vis_max_1 = m_ele, m_ele
        if (decay_type_1_i == 1):
            m_vis_min_1, m_vis_max_1 = m_muon, m_muon
        if (decay_type_1_i == 2):
            m_vis_min_1, m_vis_max_1 = m_pion, 1.5
        if (decay_type_2_i == 0):
            m_vis_min_2, m_vis_max_2 = m_ele, m_ele
        if (decay_type_2_i == 1):
            m_vis_min_2, m_vis_max_2 = m_muon, m_muon
        if (decay_type_2_i == 2):
            m_vis_min_2, m_vis_max_2 = m_pion, 1.5
        if (m_vis_1 < m_vis_min_1):
            m_vis_1 = m_vis_min_1
        if (m_vis_1 > m_vis_max_1):
            m_vis_1 = m_vis_max_1
        if (m_vis_2 < m_vis_min_2):
            m_vis_2 = m_vis_min_2
        if (m_vis_2 > m_vis_max_2):
            m_vis_2 = m_vis_max_2

        # store both tau candidate four vectors
        leg1 = ak.zip({"pt": pt_1[i], "eta": eta_1[i], "phi": phi_1[i], "mass":m_vis_1}, 
                        with_name="PtEtaPhiMLorentzVector", behavior=coffea.nanoevents.methods.nanoaod.vector.behavior)
        leg2 = ak.zip({"pt": pt_2[i], "eta": eta_2[i], "phi": phi_2[i], "mass":m_vis_2}, 
                        with_name="PtEtaPhiMLorentzVector", behavior=coffea.nanoevents.methods.nanoaod.vector.behavior)

        # store visible mass of ditau pair
        m_vis = np.sqrt(2*leg1.pt*leg2.pt*(np.cosh(leg1.eta - leg2.eta) -
                                             np.cos(leg1.phi - leg2.phi)))
        m_vis = float(m_vis[0]) #must not be an array but float or str !

        # correct initial visible masses
        if (decay_type_1_i == 2 and m_vis_1 > 1.5):
            m_vis_1 = 0.3
        if (decay_type_2_i == 2 and m_vis_2 > 1.5):
            m_vis_2 = 0.3

        # invert met covariance matrix, calculate determinant
        metcovinv_xx, metcovinv_yy = metcov_yy[i], metcov_xx[i]
        metcovinv_xy, metcovinv_yx = -metcov_xy[i], -metcov_yx[i]
        metcovinv_det = (metcovinv_xx*metcovinv_yy -
                         metcovinv_yx*metcovinv_xy)
        if (metcovinv_det<1e-10):
                print("Warning! Ill-conditioned MET covariance at event index", i)
                continue

        # perform likelihood scan
        # see http://cms.cern.ch/iCMS/jsp/openfile.jsp?tp=draft&files=AN2019_032_v3.pdf
        met_const = float(1/(2*np.pi*np.sqrt(metcovinv_det))[0])
        min_likelihood, x1_opt, x2_opt = 999, 0.01, 0.01
        mass_likelihood, met_transfer = 0, 0

        initialise = True

        # scan over weights for each ditau four-vector
        for x1 in np.arange(0.01, 1, 0.01):
            for x2 in np.arange(0.01, 1, 0.01):
                x1_min = float(min(1, np.pow((m_vis_1/m_tau), 2)))
                x2_min = float(min(1, np.pow((m_vis_2/m_tau),2)))
                if ((x1 < x1_min) or (x2 < x2_min)):
                    continue

                # test weighted four-vectors
                leg1_x1, leg2_x2 = leg1*(1/x1), leg2*(1/x2)

                ditau_test = ak.zip(
                    {
                        "px": leg1_x1.px + leg2_x2.px, 
                        "py": leg1_x1.py + leg2_x2.py, 
                        "pz": leg1_x1.pz + leg2_x2.pz, 
                        "E": leg1_x1.E + leg2_x2.E
                    }, 
                    with_name="LorentzVector") #, behavior=vector.behavior)

                nu_test = ak.zip(
                    {
                        "px": ditau_test["px"] - leg1.px - leg2.px,
                        "py": ditau_test["py"] - leg1.py - leg2.py,   
                        "pz": ditau_test["pz"] - leg1.pz - leg2.pz, 
                        "E": ditau_test["E"] - leg1.E - leg2.E  
                        
                    }
                )

                mass_ditau_test = float(np.sqrt(
                    ditau_test['E']**2 - ditau_test['px']**2 - 
                    ditau_test['py']**2 - ditau_test['pz']**2)[0])
                test_mass = mass_ditau_test

                if constrain_setting == "Window":
                    if (((test_mass < constrain_window[0]) or
                        (test_mass > constrain_window[1])) and
                        constrain):
                        continue

                # calculate mass likelihood integral
                m_shift = test_mass * delta
                if (m_shift < m_vis):
                    continue
                x1_min = float(min(1.0, np.pow((m_vis_1/m_tau),2)))
                x2_min = float(max(np.pow((m_vis_2/m_tau),2),
                             np.pow((m_vis/m_shift),2)))
                x2_max = float(min(1.0, np.pow((m_vis/m_shift),2)/x1_min))
                if (x2_max < x2_min):
                    continue
                J = float(2*np.pow(m_vis,2)) * float(np.pow(m_shift, -reg_order))
                I_x2 = float(np.log(x2_max)) - float(np.log(x2_min))
                I_tot = I_x2
                if (decay_type_1_i != 2):
                    I_m_nunu_1 = float(np.pow((m_vis/m_shift),2)) * (float(np.pow(x2_max,-1)) - float(np.pow(x2_min,-1)))
                    I_tot += I_m_nunu_1
                if (decay_type_2_i != 2):
                    I_m_nunu_2 = float(np.pow((m_vis/m_shift),2)) * I_x2 - (x2_max - x2_min)
                    I_tot += I_m_nunu_2
                mass_likelihood = 1e9 * J * I_tot

                # calculate MET transfer function
                residual_x = met_x[i] - nu_test['px']
                residual_y = met_y[i] - nu_test['pz']
                pull2 = (residual_x*(metcovinv_xx*residual_x +
                                     metcovinv_xy*residual_y) +
                         residual_y*(metcovinv_yx*residual_x +
                                     metcovinv_yy*residual_y))
                pull2_numpy = pull2.to_numpy()
                pull2_numpy /= metcovinv_det
                met_transfer = met_const*float(np.exp(-0.5*pull2_numpy))

                # calculate final likelihood, store if minimum
                likelihood = -met_transfer * mass_likelihood

                if constrain and constrain_setting == "BreitWigner":
                    mH = 125.0
                    GammaH = 0.004
                    deltaM = test_mass*test_mass - mH*mH
                    mG = test_mass*GammaH
                    BreitWigner_likelihood = 1/(deltaM*deltaM + mG*mG)
                    likelihood = likelihood*BreitWigner_likelihood

                if initialise:
                    min_likelihood = likelihood
                    x1_opt, x2_opt = x1, x2
                    initialise = False
                else:
                    if (likelihood < min_likelihood):
                        min_likelihood = likelihood
                        x1_opt, x2_opt = x1, x2

        leg1_x1, leg2_x2 = leg1*(1/x1_opt), leg2*(1/x2_opt)
        
        p4_ditau_opt = ak.zip(
                    {
                        "px": leg1_x1.px+leg2_x2.px, 
                        "py": leg1_x1.py+leg2_x2.py, 
                        "pz": leg1_x1.pz+leg2_x2.pz, 
                        "E": leg1_x1.E+leg2_x2.E
                    }, 
                    with_name="LorentzVector") #, behavior=vector.behavior)

        mass_opt = float(np.sqrt(
            p4_ditau_opt['E']**2 - p4_ditau_opt['px']**2 - 
            p4_ditau_opt['py']**2 - p4_ditau_opt['pz']**2)[0])

        pt_opt = float(np.sqrt(
            p4_ditau_opt['px']**2 + p4_ditau_opt['py']**2)[0])


        pt1_opt = float(pt_1[i]/x1_opt)
        pt2_opt = float(pt_2[i]/x2_opt)

        #from IPython import embed; embed()

        fastmttMass_values[i] = mass_opt # fastmttMass_values[i] : np.float32()
        fastmttPt_values[i] = pt_opt
        fastmttPt1_values[i] = pt1_opt
        fastmttPt2_values[i] = pt2_opt

    return fastmttMass_values, fastmttPt_values, fastmttPt1_values, fastmttPt2_values
