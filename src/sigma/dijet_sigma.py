import numpy as np
from itertools import combinations
from scipy.integrate import dblquad
from collections.abc import Callable

from setup.process_vars import SQRT_S, S, FLAVOURS, QUARKS, ANTIQUARKS, ISOSPIN_MAP, Y_MAX, PT_MAX, PT_MIN, A, Z
from setup.partonic_sigma import *
from setup.load_pdf import pdf_proton, pdf_nucleus
from setup.alpha_s import alpha_s


def dijet_sigma_dy1dy2dpt2(pt:float, y1:float, y2:float, type:str="pp", nucleon:str="p") -> tuple:
    """
    Function calculating the single parton scattering contribution to the process pA-> dijet, depending on the jet transverse momentum and the corresponding rapidities.
    In the process of the calculation, a dictionary is created with labels corresponding to the specific subprocesses being considered. This way, this function may be used
    in the calculation of the DPS contribution, since a proper combinatorial factor needs to be taken into account.

    :param pt: jet transverse momentum
    :param y1: rapidity of the first jet in the pair 
    :param y2: rapidity of the second jet in the pair
    :param type: string representing the type of the interaction, i.e. either "pp" (default) or "pPb"
    :param nucleon: character ("p" (default) or "n") representing the nucleon participating in the proton-nucleus SPS interaction; using this convention, the function calculating the DPS can reuse this functionality for each of the terms in the final expression

    :returns results: dictionary of SPS differential cross section of pN -> dijet with respect to dy1 dy2 dpt^2 in GeV^{-4}
    :returns results_sum: total SPS cross section used for calculations in GeV
    """
        
    # Calculate the momenta fractions
    x1 = pt/SQRT_S * (np.exp(y1) + np.exp(y2))
    x2 = pt/SQRT_S * (np.exp(-y1) + np.exp(-y2))

    # Handle x limits; if unphysical return an empty dict and null value
    if x1 >= 1 or x2 >= 1 or x1 <= 0 or x2 <= 0: return {}, 0.0

    # Calculate the partonic Mandelstam variables
    shat = x1 * x2 * S
    that = -pt * pt * (1 + np.exp(-(y1 - y2)))
    uhat = -pt * pt * (1 + np.exp(y1 - y2))

    # Set the hard scale by assuming the transverse momentum of the jet 
    Q2 = pt**2
    Q = pt

    # Load PDFs as dictionaries indexed with MC names
    f1_proton = {p: pdf_proton(p, x1, Q) for p in FLAVOURS}
    f2_proton = {p: pdf_proton(p, x2, Q) for p in FLAVOURS}
    f2_nucleus = {p: pdf_nucleus(p, x2, Q) for p in FLAVOURS}

    # Initialise the results dict
    results = {}

    # Calculate the cross section pre-factor 
    prefactor = np.pi * alpha_s(Q2)**2/shat**2

    # Define an incrementing function to add individual contributions depending on the type of collision
    def add(key:tuple, a:str, b:str, sigma:Callable):
        # If identical initial state partons
        if a == b:
            if type == "pp": val = x1 * f1_proton[a] * x2 * f2_proton[b] * sigma(shat, that, uhat)
            elif type == "pPb": 
                if nucleon == "p": val = x1 * f1_proton[a] * x2 * f2_nucleus[b] * sigma(shat, that, uhat)
                elif nucleon == "n": val = x1 * f1_proton[a] * x2 * f2_nucleus[ISOSPIN_MAP[b]] * sigma(shat, that, uhat)
        # else, sum both orderings but use a unique dict key to allow for proper comparison
        else: 
            if type == "pp": val = (x1 * f1_proton[a] * x2 * f2_proton[b] + x1 * f1_proton[b] * x2 * f2_proton[a]) * sigma(shat, that, uhat)
            elif type == "pPb": 
                if nucleon == "p": val = (x1 * f1_proton[a] * x2 * f2_nucleus[b] + x1 * f1_proton[b] * x2 * f2_nucleus[a]) * sigma(shat, that, uhat)
                elif nucleon == "n": val = (x1 * f1_proton[a] * x2 * f2_nucleus[ISOSPIN_MAP[b]] + x1 * f1_proton[b] * x2 * f2_nucleus[ISOSPIN_MAP[a]]) * sigma(shat, that, uhat)

        # Add the result to the dict with default value 0.0 if entry not yet present
        results[key] = results.get(key, 0.0) + val


    # --- Add the individual processes ----
    # gg -> gg (1)
    add(("gg", "gg"), "g", "g", siggggg)

    # gg -> qqbar (3)
    for q, qbar in zip(QUARKS, ANTIQUARKS): add(("gg", (q,qbar)), "g", "g", sigggqibi)

    # qqbar -> gg (3)
    for q, qbar in zip(QUARKS, ANTIQUARKS): add(((q,qbar), "gg"), q, qbar, sigqibigg)

    # gq -> gq (6)
    for q in QUARKS + ANTIQUARKS: add((("g", q), ("g", q)), "g", q, siggqigqi)

    # qq -> qq (incl. qbar qbar -> qbar qbar) (t- and u-channels) (6)
    for q in QUARKS + ANTIQUARKS: add(((q,q), (q,q)), q, q, sigqiqiqiqi)

    # --- q q' -> q q' (t-channel) (12) ---
    # q q' -> q q' for quarks (3)
    for q1, q2 in combinations(QUARKS, 2): add(((q1,q2), (q1,q2)), q1, q2, sigqiqjqiqj)

    # q q' -> q q' for antiquarks (3)
    for q1, q2 in combinations(ANTIQUARKS, 2): add(((q1,q2), (q1,q2)), q1, q2, sigqiqjqiqj)
                
    # q q' -> q q' for quark and antiquark (6)
    for q in QUARKS:
        for qbar in ANTIQUARKS:
            if q + "bar" == qbar: continue # skip if same flavour
            add(((q,qbar), (q,qbar)), q, qbar, sigqiqjqiqj)
    # -------------------

    # q qbar -> q' qbar' (6)
    for q, qbar in zip(QUARKS, ANTIQUARKS):
        for q2, q2bar in zip(QUARKS, ANTIQUARKS):
            if q2 == q: continue # skip if same flavour
            add(((q,qbar), (q2,q2bar)), q, qbar, sigqibiqjbj)

    # q qbar -> q qbar (s- + t-channels) (3)
    for q, qbar in zip(QUARKS, ANTIQUARKS): add(((q,qbar), (q,qbar)), q, qbar, sigqibiqibi)
    # ------------------------------------------------

    # Rescale the cross sections by the prefactor
    results.update((key, value * prefactor) for key, value in results.items())

    # Sum the cross section contributions for SPS
    results_sum = sum(results.values())

    return results, results_sum
    

def dijet_sigma_dDelta_y(Delta_y:float, type:str="pp", epsabs:float=1e-4, epsrel:float=1e-4) -> float:
    """
    Integrate the differential dijet production cross section over Y and p_T^2 to obtain the expression in terms of jet rapidity differences.
    Introduces rapidity difference Delta_y = y_1 - y_2 and an extra auxiliary variable Y = 1/2 * (y_1 + y_2).

                    dsigma/dDelta_y = int dp_T^2 int dY dsigma/dy_1dy_2dp_T^2 = int dp_T 2 * p_T int dY dsigma/dy_1dy_2dp_T^2
    
    The limits for p_T^2 integration are entirely determined on PT_MIN and PT_MAX.
    The limits for Y integration are set through the assumption that |y_1| < y_max and |y_2| < y_max, from where the intersection implies
                    -(y_max - |Delta y|/2) < Y < (y_max - |Delta y|/2)

    :param Delta_y: rapidity difference
    :param type: string representing the type of the interaction, i.e. either "pp" (default) or "pPb"
    :param epsabs: absolute integration precision for dblquad (default = 1e-4)
    :param epsrel: relative integration precision for dblquad (default = 1e-4)

    :returns integral[0]: result of the double integration 
    """

    # Set limits for Y integration
    Y_min = -Y_MAX + np.abs(Delta_y)/2
    Y_max = Y_MAX - np.abs(Delta_y)/2

    # Set limits for pT integration
    pt_min = PT_MIN
    pt_max = PT_MAX

    # Define the inner integrand function
    def integrand(Y:float, pt:float) -> float:
        y1 = Y + Delta_y/2
        y2 = Y - Delta_y/2

        # Get the sum over dictionary entries
        if type == "pp": return 2 * pt * dijet_sigma_dy1dy2dpt2(pt, y1, y2, type="pp", nucleon="p")[1]
        # Z protons and (A-Z) neutrons contribution from the nucleus; use the isospin map
        elif type == "pPb": return 2 * pt * (Z * dijet_sigma_dy1dy2dpt2(pt, y1, y2, type="pPb", nucleon="p")[1] + (A-Z) * dijet_sigma_dy1dy2dpt2(pt, y1, y2, type="pPb", nucleon="n")[1])
        else: raise ValueError("Unknown collision type parameter...")

    # Integrate outer integral over pt
    integral = dblquad(integrand, pt_min, pt_max, Y_min, Y_max, epsabs=epsabs, epsrel=epsrel)

    # Return the result dsigma/dDelta_y in GeV^{-2}
    return integral[0]
