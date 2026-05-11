from double_dijet_sigma import double_dijet_sigma_dy1dy2dy3dy4pt12pt22, nuclear_double_dijet_sigma_dy1dy2dy3dy4pt12pt22
from process_vars import PT_MIN, PT_MAX, Y_MAX, RADIUS, CONV_GEV_NB, A

import numpy as np


def injet_double_dijet_sigma_total(N:int, type:str):
    """
    Evaluate the cross section for jet 3 to lie within the cone of radius R defined by jet 1, taking advantage of the Iterated Cone with Progressive Removal of Particles (ICPR) algo definition of the jet cone through the following formula
                
                Delta R_is^2 = (y_i - y_s)^2 + (phi_i - phi_s)^2, where y_s is the rapidity of the source jet, and phi_s is its azimuthal angle coordinate.
    This can be effectively enforced by using a Heaviside function within the integrand. For that purpose, the cross section is averaged over the azimuths to obtain a differential in phi_3 and phi_1 
                1/(2pi)^2 int dphi_3 dphi_1 dsigma/dy1dy3 Theta(R^2 - delta_y^2 - delta_phi^2).

    From there for the integration over the two rapidities and two angles, a change of variables may be performed
                Y = 1/2 (y1 + y3), delta_y = y1 - y3, Phi = 1/2 (phi_1 + phi_3), delta_phi = phi_1 - phi_3

    The full differential cross section as obtained by double_dijet_sigma_dy1dy2dy3dy4dpt12dpt22 is integrated, using a Monte Carlo method, over the unconstrained variables pt12, pt12, y2, and y4 and the relevant constrained variables Y, Phi, delta_y, delta_phi.
    After the full integration, the total cross will be modulated by an extra factor due to the geometrical constaints on the jet 1 + jet 3 kinematics.

    :param N: number of random points to sample for each parameter
    :param type: string representing the type of the interaction, i.e. either "pp" (default) or "pPb"

    :returns result: Monte Carlo value of the cross section in nb
    :returns err: error estimate of the cross section in nb
    """

    # Sample the phase space in a random fashion
    pt1 = np.random.uniform(PT_MIN, PT_MAX, N)
    pt2 = np.random.uniform(PT_MIN, PT_MAX, N)

    y2 = np.random.uniform(-Y_MAX, Y_MAX, N)
    y4 = np.random.uniform(-Y_MAX, Y_MAX, N)
    
    Delta_y = np.random.uniform(-RADIUS, RADIUS, N)
    
    # Note a uniform range -- dependent on Delta_y
    Y_lim = Y_MAX - np.abs(Delta_y)/2
    Y = np.random.uniform(-1, 1, N) * Y_lim

    # Calculate the phase space volume of integration: over pt12, pt22, y2, y4, Delta_y
    volume = (PT_MAX - PT_MIN)**2 * (2 * Y_MAX)**2 * (2 * RADIUS)

    # Calculate the extra geometrical factor due to Heaviside constraints in rapidity and angular integration
    geo = 1/np.pi * np.sqrt(RADIUS**2 - Delta_y**2) - (RADIUS**2 - Delta_y**2) / (4 * np.pi**2)

    # Initialise the arrays of the integrand results
    dsigma = np.zeros(N)

    # Change variables from integration over dy1dy3
    y1 = Y + Delta_y/2
    y3 = Y - Delta_y/2

    # Use pp or pPb functions to evaluate the integrand, depending on the type parameter
    if type == "pp":
        for i in range(N): dsigma[i] = 2 * pt1[i] * 2 * pt2[i] * double_dijet_sigma_dy1dy2dy3dy4pt12pt22(pt1[i], pt2[i], y1[i], y2[i], y3[i], y4[i], type="pp", nucleon_1="p", nucleon_2="p")
    elif type == "pPb":
        for i in range(N): dsigma[i] = 2 * pt1[i] * 2 * pt2[i] * nuclear_double_dijet_sigma_dy1dy2dy3dy4pt12pt22(pt1[i], pt2[i], y1[i], y2[i], y3[i], y4[i])
    else: raise ValueError("Unknown collision type parameter...")

    # Include geometry and the Delta_y-dependent volume factor for integration over Y
    dsigma = dsigma * geo * 2 * Y_lim

    # Assign a compensation factor due to jet permutations 1 <-> 2; this ensures that no double-counting is present in the Monte Carlo phase space
    comp = 1/2

    # Estimate the Monte Carlo cross section expressed in nb
    result = volume * np.mean(dsigma, dtype=float) * CONV_GEV_NB * comp
    
    # Estimate the error in the result
    err = volume * comp / np.sqrt(N) * np.std(dsigma) * CONV_GEV_NB

    return result, err