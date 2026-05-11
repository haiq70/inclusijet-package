from dijet_sigma import dijet_sigma_dy1dy2dpt2
from partonic_sigma import *
from process_vars import SIGMA_EFF, PT_MIN, PT_MAX, Y_MAX, CONV_GEV_NB, A, Z, RA

import numpy as np

def double_dijet_sigma_dy1dy2dy3dy4pt12pt22(pt1:float, pt2:float, y1:float, y2:float, y3:float, y4:float, type:str="pp", nucleon_1:str="p", nucleon_2:str="p") -> float:
    """
    Compute the influence of double parton scattering to pA -> 4 jets using the factorisation ansatz. The implementation takes advantage of the dictionaries built in the calculation of the SPS cross section, to compare the 
    individual processes via their labels and assign a combinatorial factor of 1/2 to the total cross section if both the initial and final states are the same. Else, the factor remains unity.

    :param pt1: transverse momentum of the first pair of jets
    :param pt2: transverse momentum of the second pair of jets
    :param y1: jet 1 rapidity
    :param y2: jet 2 rapidity
    :param y3: jet 3 rapidity
    :param y4: jet 4 rapidity
    :param type: string representing the type of the interaction, i.e. either "pp" or "pPb"; default setting corresponds to proton-proton collisions
    :param nucleon_1: character ("p" or "n") representing the first nucleon participating in the proton-nucleus SPS interaction; using this convention, the function calculating the DPS can reuse this functionality for each of the terms in the final expression
    :param nucleon_2: character ("p" or "n") representing the second nucleon

    :returns total (/ SIGMA_EFF): DPS differential cross section of pp -> 4 jets with respect to dy1 dy2 dy3 dy4 dpt_1^2 dpt_2^2 in GeV^{-6}
    """

    # Extract the SPS dicts
    sub1, _ = dijet_sigma_dy1dy2dpt2(pt1, y1, y2, type, nucleon_1)
    sub2, _ = dijet_sigma_dy1dy2dpt2(pt2, y3, y4, type, nucleon_2)
    
    
    # Loop through each subprocess
    total = 0.0
    for key1, val1 in sub1.items():
        for key2, val2 in sub2.items():
            # Apply the combinatorial factor: 1/2 if both initial and final states are identical
            if key1 == key2: C = 0.5
            else: C = 1.0

            total += C * val1 * val2

    if type == "pp": return total / SIGMA_EFF
    # If proton-nucleus collision indicated, return the non-scaled DPS contribution, corresponding to the product of two pA -> 2 jets
    elif type == "pPb": return total
    
    else: raise ValueError("Unknown collision type parameter...")

def nuclear_double_dijet_sigma_dy1dy2dy3dy4pt12pt22(pt1:float, pt2:float, y1:float, y2:float, y3:float, y4:float) -> float:
    """
    Evaluate the differential cross section for the double parton scattering pA -> 4 jets using modelling in arXiv:2001.04256v2, assuming no spatial dependence of the nuclear PDFs and using the hard-sphere model for the nucleon density.
    This implementation takes advantage of the above function double_dijet_sigma_dy1dy2dy3dy4pt12pt22 to compute the products of SPS scatterings.

    :param pt1: transverse momentum of the first pair of jets
    :param pt2: transverse momentum of the second pair of jets
    :param y1: jet 1 rapidity
    :param y2: jet 2 rapidity
    :param y3: jet 3 rapidity
    :param y4: jet 4 rapidity

    :returns: sum of nucleon contribution to the cross section, representing the approximation of the DPS effect in proton-lead collisions measured in GeV^{-6} (due to differentials in transverse momentum)
    """

    # Compute the constant nucleus-dependent factor
    A_factor = (A - 1)/A * 9/(8 * np.pi * RA**2)
    
    # Extract the subsequent nucleon contributions
    # Two partons from one single nucleon inside the nucleus -- 2x proton or 2x neutron
    sigma_same_pp = double_dijet_sigma_dy1dy2dy3dy4pt12pt22(pt1, pt2, y1, y2, y3, y4, type="pPb", nucleon_1="p", nucleon_2="p") * Z
    sigma_same_nn = double_dijet_sigma_dy1dy2dy3dy4pt12pt22(pt1, pt2, y1, y2, y3, y4, type="pPb", nucleon_1="n", nucleon_2="n") * (A-Z)

    # Two partons from two different nucleons inside the nucleus -- 2x proton, 2x neutron, proton-neutron, neutron-proton
    sigma_diff_pp = double_dijet_sigma_dy1dy2dy3dy4pt12pt22(pt1, pt2, y1, y2, y3, y4, type="pPb", nucleon_1="p", nucleon_2="p") * Z * (Z-1)
    sigma_diff_nn = double_dijet_sigma_dy1dy2dy3dy4pt12pt22(pt1, pt2, y1, y2, y3, y4, type="pPb", nucleon_1="n", nucleon_2="n") * (A-Z) * (A-Z-1)
    # Note that for identical final states, the cross terms are the same
    sigma_diff_np = double_dijet_sigma_dy1dy2dy3dy4pt12pt22(pt1, pt2, y1, y2, y3, y4, type="pPb", nucleon_1="n", nucleon_2="p") * Z * (A-Z)
    sigma_diff_pn = double_dijet_sigma_dy1dy2dy3dy4pt12pt22(pt1, pt2, y1, y2, y3, y4, type="pPb", nucleon_1="p", nucleon_2="n") * Z * (A-Z)
    

    # Return the total cross section in GeV^{-6}
    return (sigma_same_pp + sigma_same_nn) * (1/SIGMA_EFF + A_factor) + (sigma_diff_pp + sigma_diff_nn + sigma_diff_pn + sigma_diff_np) * A_factor


def double_dijet_sigma_Delta_y_max(N:int, n_bins:int, y_range_min:float, y_range_max:float, type:str="pp") -> tuple:
    """
    Monte Carlo integration of the differential cross section calculated by double_dijet_sigma_dy1dy2dy3dy4pt12pt22 to express the differential in terms of the rapidity distance between the two jets most separated in rapidity space.
    The function first samples random values for the two transverse momenta and rapidities, based on the limits set by the global variables. Next, for each of these samples, it determines the maximum rapidity distance between jets and computes the integrand arising when changing coordinates
    to ones expressing the maximum separation. Using numpy's functionalities, it then builds a histogram based on the Delta_y_max and dsigma data. The return value is the Monte Carlo estimate of the integral, based on the available phase space. 

    :param N: number of random points to sample for each parameter
    :param n_bins: number of histogram bins to produce
    :param y_range_min: minimum value for the Delta_y range, used for histogram building
    :param y_range_max: maximum value for the Delta_y range, used for histogram building
    :param type: string representing the type of the interaction, i.e. either "pp" (default) or "pPb"

    :returns comp * result * CONV_GEV_NB: binned values of the cross-section (y-data) in nb, compensated for the identical phase space factors 
    :returns bin_centres: bin centre (x-data) locations 
    """
    
    # Monte Carlo sampling based on the global limits
    pt1 = np.random.uniform(PT_MIN, PT_MAX, N)
    pt2 = np.random.uniform(PT_MIN, PT_MAX, N)

    y1 = np.random.uniform(-Y_MAX, Y_MAX, N)
    y2 = np.random.uniform(-Y_MAX, Y_MAX, N)
    y3 = np.random.uniform(-Y_MAX, Y_MAX, N)
    y4 = np.random.uniform(-Y_MAX, Y_MAX, N)

    # Initialise the x and y values lists
    Delta_y_max_vals, dsigma_vals = [], []
    
    # Iterate over each sampled point
    for i in range(N):
        # Calculate the largest rapidity separation between all possible combinations
        Delta_y_max = max(abs(y1[i]-y2[i]), abs(y1[i]-y3[i]), abs(y1[i]-y4[i]), abs(y2[i]-y3[i]), abs(y2[i]-y4[i]), abs(y3[i]-y4[i]))
        Delta_y_max_vals.append(Delta_y_max)

        # Calculate the integrand
        if type == "pp": dsigma = 2 * pt1[i] * 2 * pt2[i] * double_dijet_sigma_dy1dy2dy3dy4pt12pt22(pt1[i], pt2[i], y1[i], y2[i], y3[i], y4[i], type="pp", nucleon_1="p", nucleon_2="p")
        elif type == "pPb": dsigma = 2 * pt1[i] * 2 * pt2[i] * nuclear_double_dijet_sigma_dy1dy2dy3dy4pt12pt22(pt1[i], pt2[i], y1[i], y2[i], y3[i], y4[i])
        else: raise ValueError("Unknown collision type parameter...")

        dsigma_vals.append(dsigma)

    # Convert to numpy arrays
    Delta_y_max_vals = np.array(Delta_y_max_vals)
    dsigma_vals = np.array(dsigma_vals)

    # Build a histogram in the specified range of max rapidity differences 
    bins = np.linspace(y_range_min, y_range_max, n_bins+1)
    hist, edges = np.histogram(Delta_y_max_vals, bins=bins, weights=dsigma_vals, density=False)

    bin_width = edges[1] - edges[0]
    bin_centres = 1/2 * (bins[:-1] + bins[1:])

    # Calculate the phase space volume due to the MC integration scheme; two (PT_MIN, PT_MAX) ranges and four (-Y_MAX, Y_MAX) ranges 
    volume = (PT_MAX - PT_MIN)**2 * (2 * Y_MAX)**4

    # Calculate the value of the integral evaluated at each histogram Delta_y_max value
    result = volume * hist / (N * bin_width)

    # Assign a compensation factor due to jet permutations 1 <-> 2; this ensures that no double-counting is present in the Monte Carlo phase space
    comp = 1/2

    # Return the compensated cross section in nb and bin centre locations
    return comp * result * CONV_GEV_NB, bin_centres


def double_dijet_sigma_total(N:int, type:str="pp") -> tuple:
    """
    Performing a complete Monte Carlo integration over the rapidity and transverse momentum differentials in dsigma/dy1dy2dy3dy4dpt12dpt22, returning the total cross section for the DPS of two protons or a proton and a neutron into four jets.
    The algorithm is precisely the same as in double_dijet_sigma_Delta_y_max, but instead of producing a histogram in terms of rapidity differences, the function simply calculates the mean of the integrand samples, normalised through the number of points and the integration phase space. 
    Note that the range of the integration is set by the global parameters, defined in the process_vars.py file. 

    :param N: number of random points to sample for each parameter
    :param type: string representing the type of the interaction, i.e. either "pp" (default) or "pPb"

    :returns result: Monte Carlo value of the cross section in nb
    :returns err: error estimate of the cross section in nb
    """

    # Sample the phase space in a random fashion
    pt1 = np.random.uniform(PT_MIN, PT_MAX, N)
    pt2 = np.random.uniform(PT_MIN, PT_MAX, N)

    y1 = np.random.uniform(-Y_MAX, Y_MAX, N)
    y2 = np.random.uniform(-Y_MAX, Y_MAX, N)
    y3 = np.random.uniform(-Y_MAX, Y_MAX, N)
    y4 = np.random.uniform(-Y_MAX, Y_MAX, N)

    # Calculate the phase space volume of integration
    volume = (PT_MAX - PT_MIN)**2 * (2 * Y_MAX)**4 

    # Initialise the arrays of the integrand results
    dsigma = np.zeros(N)

    # Use pp or pPb functions to evaluate the integrand, depending on the type parameter
    if type == "pp":
        for i in range(N): dsigma[i] = 2 * pt1[i] * 2 * pt2[i] * double_dijet_sigma_dy1dy2dy3dy4pt12pt22(pt1[i], pt2[i], y1[i], y2[i], y3[i], y4[i], type="pp", nucleon_1="p", nucleon_2="p")
    elif type == "pPb":
        for i in range(N): dsigma[i] = 2 * pt1[i] * 2 * pt2[i] * nuclear_double_dijet_sigma_dy1dy2dy3dy4pt12pt22(pt1[i], pt2[i], y1[i], y2[i], y3[i], y4[i])
    else: raise ValueError("Unknown collision type parameter...")

    # Assign a compensation factor due to jet permutations 1 <-> 2; this ensures that no double-counting is present in the Monte Carlo phase space
    comp = 1/2

    # Estimate the Monte Carlo cross section expressed in nb
    result = volume * comp *  np.mean(dsigma, dtype=float) * CONV_GEV_NB

    # Estimate the error in the result
    err = volume * comp / np.sqrt(N) * np.std(dsigma) * CONV_GEV_NB

    return result, err