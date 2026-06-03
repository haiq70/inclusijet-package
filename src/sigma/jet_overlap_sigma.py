from sigma.double_dijet_sigma import double_dijet_sigma_dy1dy2dy3dy4pt12pt22, nuclear_double_dijet_sigma_dy1dy2dy3dy4pt12pt22
from setup.process_vars import PT_MIN, PT_MAX, Y_MAX, RADIUS, CONV_GEV_NB, A

import numpy as np

def injet_double_dijet_sigma_dPT(PT:float, N:int, type:str, pt_cut:float, radius:float) -> tuple:
    """
    Evaluating the differential (w.r.t total reconstructed jet transverse momentum) cross section of double parton scattering production of four-jets, represented in terms of the azimuthal angle and transverse momentum of the reconstructed jet, restricted to a condition where one of the jets lies within the cone developed by another.
    This kinematical constraint is enforced by an introduction of a Heaviside function, relating the radius of the cone, rapidity separation of the two relevant jets, as well as their azimuthal distance. The full transformation resulting in the final form of this cross section may be found in the project report.
    The approach to multidimensional integration of the cross section follows the proprietary Monte Carlo scheme, developed in the previous cross section implementations.

    :param PT: total jet transverse momentum as determined by the experimental measurement
    :param N: number of random points to sample for each parameter
    :param type: string representing the type of the interaction, i.e. either "pp" (default) or "pPb"
    :param pt_cut: minimum transverse momentum of a jet pair, used as a cut in the cross section evaluation
    :param radius: parametrically small radius of the jet cone, as defined by jet labelled by the index 1
    
    :returns result: Monte Carlo value of the cross section in nb
    :returns err: error estimate of the cross section in nb 
    """
    
    # Sample the integration variables    
    # (Unconstrained) Rapidities 
    y2 = np.random.uniform(-Y_MAX, Y_MAX, N)
    y4 = np.random.uniform(-Y_MAX, Y_MAX, N)
    
    # (Unconstrained) Auxiliary momentum of the reconstructed jet
    qt = np.random.uniform(0, PT_MAX, N)

    # Rapidity separation between two relevant jets (to be constrained using Heaviside)
    Delta_y = np.random.uniform(-2*Y_MAX, 2*Y_MAX, N)
    
    # (Unconstrained) Auxiliary rapidity variable with limits dependent on Delta_y
    Y_lim = Y_MAX - np.abs(Delta_y)/2
    Y = np.random.uniform(-1, 1, N) * Y_lim

    # Get the constrained rapidities after change of variables
    y1 = Y + Delta_y/2
    y3 = Y - Delta_y/2

    # Azimuthal angle of the auxiliary transverse momentum variable (to be constrained using Heaviside)
    phi = np.random.uniform(0, 2*np.pi, N)

    # Calculate the volume of integration phase space, excluding the terms dependent on the integration variables' values
    volume = (2 * Y_MAX)**2 * (PT_MAX) * (4*Y_MAX) * (2*np.pi)

    # Evaluate the integrands for the specified type of interactions
    dsigma = np.zeros(N, dtype=float)

    # Calculate the angles and momenta of the jet pairs, determined by momentum conservation condition
    ptA = np.sqrt(0.25 * PT**2 + qt**2 + PT * qt * np.cos(phi))
    ptB = np.sqrt(0.25 * PT**2 + qt**2 - PT * qt * np.cos(phi))

    # Ensure proper treatment of branch cuts in the inverse tangent
    phi1 = np.arctan2(2 * qt * np.sin(phi), PT + 2 * qt * np.cos(phi)) # (-pi, pi)
    phi3 = np.arctan2((-1) * 2 * qt * np.sin(phi), PT - 2 * qt * np.cos(phi)) # (-pi, pi)
    # Wrap the separation to include values in the range (-pi, pi); correctly handle jets close to the branch cut
    Delta_phi = (phi1 - phi3 + np.pi) % (2*np.pi) - np.pi

    # Cast Heaviside condition as a mask on the numpy arrays
    heaviside = (radius**2 - Delta_y**2 - Delta_phi**2) > 0
    # Include a physical cut on the momenta to eliminate singular behaviour of the cross section and account for detector acceptance 
    physical = (ptA >= pt_cut) & (ptB >= pt_cut) & (ptA <= PT) & (ptB <= PT)

    # Get index of array for iteration
    idx = np.where(heaviside & physical)[0]
    print(len(idx)) # Debug printing

    # Calculate only the values of the integrand which are nonzero; attempt to make the iteration more computationally efficient 
    for i in idx:
        # Calculate the integrand, with pre-factors and Delta_y-dependent volume factor included  
        if type == "pp": dsigma[i] = 2/(np.pi) * qt[i] * (2*Y_lim[i]) * double_dijet_sigma_dy1dy2dy3dy4pt12pt22(float(ptA[i]), float(ptB[i]), float(y1[i]), float(y2[i]), float(y3[i]), float(y4[i]), type="pp", nucleon_1="p", nucleon_2="p")
        elif type == "pPb": dsigma[i] = 2/(np.pi) * qt[i] * (2*Y_lim[i]) * nuclear_double_dijet_sigma_dy1dy2dy3dy4pt12pt22(float(ptA[i]), float(ptB[i]), float(y1[i]), float(y2[i]), float(y3[i]), float(y4[i]))
        else: raise ValueError("Unknown interaction type parameter...")

    # Assign a compensation factor due to jet permutations 1 <-> 2; this ensures that no double-counting is present in the Monte Carlo phase space
    comp = 1/2

    # Estimate the Monte Carlo cross section expressed in GeV^{-4} (1/PT dsigma/dPT)
    result = volume * np.mean(dsigma, dtype=float) * comp
    
    # Estimate the error in the result
    err = volume * comp / np.sqrt(N) * np.std(dsigma)

    return result, err


def injet_double_dijet_sigma_total(N:int, type:str, pt_cut:float, radius:float) -> tuple:
    """
    Evaluating the total cross section of double parton scattering production of four-jets, represented in terms of the azimuthal angle and transverse momentum of the reconstructed jet, restricted to a condition where one of the jets lies within the cone developed by another.
    This kinematical constraint is enforced by an introduction of a Heaviside function, relating the radius of the cone, rapidity separation of the two relevant jets, as well as their azimuthal distance. The full transformation resulting in the final form of this cross section may be found in the project report.
    The approach to multidimensional integration of the cross section follows the proprietary Monte Carlo scheme, developed in the previous cross section implementations.

    :param N: number of random points to sample for each parameter
    :param type: string representing the type of the interaction, i.e. either "pp" (default) or "pPb"
    :param pt_cut: minimum transverse momentum of a jet pair, used as a cut in the cross section evaluation
    :param radius: parametrically small radius of the jet cone, as defined by jet labelled by the index 1
    
    :returns result: Monte Carlo value of the cross section in nb
    :returns err: error estimate of the cross section in nb 
    """
    
    # Sample the integration variables
    # Observable jet variables
    PT = np.random.uniform(PT_MIN, PT_MAX, N)
    
    # (Unconstrained) Rapidities 
    y2 = np.random.uniform(-Y_MAX, Y_MAX, N)
    y4 = np.random.uniform(-Y_MAX, Y_MAX, N)
    
    # (Unconstrained) Auxiliary momentum of the reconstructed jet
    qt = np.random.uniform(0, PT_MAX, N)

    # Azimuthal angle of the auxiliary transverse momentum variable (to be constrained using Heaviside)
    phi = np.random.uniform(0, 2*np.pi, N)

    # Rapidity separation between two relevant jets (to be constrained using Heaviside)
    Delta_y = np.random.uniform(-2*Y_MAX, 2*Y_MAX, N)
    
    # (Unconstrained) Auxiliary rapidity variable with limits dependent on Delta_y
    Y_lim = Y_MAX - np.abs(Delta_y)/2
    Y = np.random.uniform(-1, 1, N) * Y_lim

    # Get the constrained rapidities after change of variables
    y1 = Y + Delta_y/2
    y3 = Y - Delta_y/2

    # Calculate the volume of integration phase space, excluding the terms dependent on the integration variables' values
    volume = (2 * Y_MAX)**2 * (4*Y_MAX) * (2*np.pi) * (PT_MAX - PT_MIN) * (PT_MAX)

    # Evaluate the integrands for the specified type of interactions
    dsigma = np.zeros(N, dtype=float)

    # Calculate the angles and momenta of the jet pairs, determined by momentum conservation condition
    ptA = np.sqrt(0.25 * PT**2 + qt**2 + PT * qt * np.cos(phi))
    ptB = np.sqrt(0.25 * PT**2 + qt**2 - PT * qt * np.cos(phi))

    # Ensure proper treatment of branch cuts in the inverse tangent
    phi1 = np.arctan2(2 * qt * np.sin(phi), PT + 2 * qt * np.cos(phi)) # (-pi, pi)
    phi3 = np.arctan2((-1) * 2 * qt * np.sin(phi), PT - 2 * qt * np.cos(phi)) # (-pi, pi)
    # Wrap the separation to include values in the range (-pi, pi); correctly handle jets close to the branch cut
    Delta_phi = (phi1 - phi3 + np.pi) % (2*np.pi) - np.pi 

    # Cast Heaviside condition as a mask on the numpy arrays
    heaviside = (radius**2 - Delta_y**2 - Delta_phi**2) > 0
    # Include a physical cut on the momenta to eliminate singular behaviour of the cross section and account for detector acceptance 
    physical = (ptA >= pt_cut) & (ptB >= pt_cut) & (ptA <= PT) & (ptB <= PT)

    # Get index of array for iteration
    idx = np.where(heaviside & physical)[0]
    print(len(idx))

    # Calculate only the values of the integrand which are nonzero; attempt to make the iteration more computationally efficient 
    for i in idx:
        # Calculate the integrand, with pre-factors and Delta_y-dependent volume factor included  
        if type == "pp": dsigma[i] = 2/(np.pi) * PT[i] * qt[i] * (2*Y_lim[i]) * double_dijet_sigma_dy1dy2dy3dy4pt12pt22(float(ptA[i]), float(ptB[i]), float(y1[i]), float(y2[i]), float(y3[i]), float(y4[i]), type="pp", nucleon_1="p", nucleon_2="p")
        elif type == "pPb": dsigma[i] = 2/(np.pi) * PT[i] * qt[i] * (2*Y_lim[i]) * nuclear_double_dijet_sigma_dy1dy2dy3dy4pt12pt22(float(ptA[i]), float(ptB[i]), float(y1[i]), float(y2[i]), float(y3[i]), float(y4[i]))
        else: raise ValueError("Unknown interaction type parameter...")

    # Assign a compensation factor due to jet permutations 1 <-> 2; this ensures that no double-counting is present in the Monte Carlo phase space
    comp = 1/2

    # Estimate the Monte Carlo cross section expressed in nb
    result = volume * np.mean(dsigma, dtype=float) * CONV_GEV_NB * comp
        
    # Estimate the error in the result
    err = volume * comp / np.sqrt(N) * np.std(dsigma) * CONV_GEV_NB

    return result, err