import numpy as np
import matplotlib.pyplot as plt
import logging
from scipy.optimize import curve_fit
from multiprocessing import Pool, cpu_count

from sigma.dijet_sigma import dijet_sigma_dDelta_y, dijet_sigma_dpt, dijet_sigma_total
from sigma.double_dijet_sigma import double_dijet_sigma_Delta_y_max, double_dijet_sigma_total
from sigma.jet_overlap_sigma import injet_double_dijet_sigma_total, injet_double_dijet_sigma_dPT
from setup.process_vars import Y_MAX, PT_MIN, PT_MAX, SQRT_S, PDF_PROTON_TITLE, PDF_NUCLEUS_TITLE, CONV_GEV_NB, A, RADIUS, PT_CUT

# Matplotlib config settings
plt.rcParams.update({
    "font.family":          "serif",
    "font.serif":           ["DejaVu Serif"],
    "xtick.minor.visible":  True,
    "ytick.minor.visible":  True,

})

# Set up logging info for programme execution
logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


def plot_sigma(Delta_y:np.ndarray, N:int=60000, n_bins:int=12, y_range_min:float=4.0, y_range_max:float=9.3, collision_type:str="pp") -> None:
    """
    Visualising the results of the cross section for pp -> 2 jets and DPS pp -> 4 jets for a range of rapidity differences.

    :param Delta_y: rapidity differences between two jets
    :param N: Monte Carlo samples to generate
    :param n_bins: number of histogram bins
    :param y_range_min: minimum of the range of rapidity distances to plot over (x-values)
    :param y_range_max: maximum of the range of rapidity distances to plot over (x-values)
    :param collision_type: string representing the type of the interaction, i.e. either pp ("proton") or pA ("nuclear"); default setting corresponds to proton-proton collisions
    """ 

    logger.info("Plotting SPS vs DPS sigma in pp...")

    # Run SPS
    def sps() -> np.ndarray:
        # Initialise a list of results for dsigma/dDelta_y
        dsigma = np.array([dijet_sigma_dDelta_y(i, collision_type) for i in Delta_y]) * CONV_GEV_NB
        return dsigma

    # Run DPS
    def dps() -> tuple:
        logger.info(f"Running Monte Carlo integration for DPS w/ N={N} and n_bins = {n_bins}")
        res, bin_centres = double_dijet_sigma_Delta_y_max(N, n_bins, y_range_min, y_range_max, collision_type)
        logger.info(f"DPS in pp -> 4 jets computed successfully")
        return res, bin_centres
    
    # Extract values
    dsigma = sps()
    res, bin_centres = dps() 

    # Plotting
    fig, ax = plt.subplots(figsize=(9,6))

    if type == "pp": plt.suptitle(rf"Differential cross section for ${type} \rightarrow$ jets w/{PDF_PROTON_TITLE}")
    elif type == "pPb": plt.suptitle(rf"Differential cross section for ${type} \rightarrow$ jets w/{PDF_NUCLEUS_TITLE}")
    else: raise ValueError("Unknown collision type parameter...")

    plt.title(rf"|$y_1, y_2$| < {Y_MAX}, $p_T \in [{PT_MIN}, {PT_MAX}]$ GeV" + r" @ $\sqrt{s}$ = " + f"{SQRT_S/1000} TeV")

    ax.set_yscale('log')
    ax.set_ylim(0.5*10**(0), 10**(5))

    # x-range padding
    x_eps = y_range_min/10
    ax.set_xlim(y_range_min - x_eps, y_range_max)

    ax.set_xlabel(r"Rapidity difference $\Delta y$")
    ax.set_ylabel(r"$d\sigma/d\Delta y$ [nb]")

    plt.plot(Delta_y, dsigma, label="LO pQCD dijet", linestyle='--', color='k')
    plt.scatter(bin_centres, res, color='r', label="DPS double dijets", s=16)

    # Add Monte Carlo integration setup data
    ax.text(0.03, 0.04, rf"Monte Carlo: $N={N}$ with {n_bins} bins", color="gray", transform=ax.transAxes, fontsize=9.5, verticalalignment="bottom")

    plt.legend()
    plt.grid(alpha=0.4)
    plt.show()


def plot_ratio(ratio_type:str="dps", y_range_min:float=0.0, y_range_max:float=9.0, N:int=30000, n_bins:int=10) -> None:
    """
    Define the plotting scheme for the ratio of 4 jet/2 jet production in proton-proton and proton-nucleus collisions, as a function of rapidity differences.

    :param ratio_type: string representing the type of ratio to plot ("jets" or "dps")
    :param y_range_min: minimum of the range of rapidity distances to plot over (x-values)
    :param y_range_max: maximum of the range of rapidity distances to plot over (x-values)
    :param N: Monte Carlo samples to generate
    :param n_bins: number of histogram bins
    """

    # Get cross section pp -> 4 jets
    logger.info(f"Running Monte Carlo integration for pp -> 4 jets w/ N={N}, n_bins={n_bins}...")
    int_sigma_4jets_pp, Delta_y = double_dijet_sigma_Delta_y_max(N, n_bins, y_range_min, y_range_max, type="pp")
    logger.info("pp -> 4 jets computed successfully.")

    # Get cross section pPb -> 4 jets
    logger.info(f"Runinng Monte Carlo integration for pPb -> 4 jets w/ N={N}, n_bins={n_bins}...")
    int_sigma_4jets_pPb, _ = double_dijet_sigma_Delta_y_max(N, n_bins, y_range_min, y_range_max, type="pPb")
    logger.info("pPb -> 4 jets computed successfully.")


    # Plotting
    fig, ax = plt.subplots(figsize=(9,6))
    plt.title(rf"|$y_1, y_2$| < {Y_MAX}, $p_T \in [{PT_MIN}, {PT_MAX}]$ GeV" + r" @ $\sqrt{s}$ = " + f"{SQRT_S/1000} TeV")
    ax.set_xlabel(r'Rapidity difference $\Delta y$')

    # if 4-jet/2-jet production to be compared
    if ratio_type == "jets":
        logger.info("Plotting 4-jet/2-jet ratio in pp and pPb")
        # Get cross section pp -> 2 jets in nb for the already determined Delta_y bin centres
        int_sigma_2jets_pp = np.array([dijet_sigma_dDelta_y(dy, type="pp") for dy in Delta_y]) * CONV_GEV_NB

        # Get cross section pPb -> 2 jets
        int_sigma_2jets_pPb = np.array([dijet_sigma_dDelta_y(dy, type="pPb") for dy in Delta_y]) * CONV_GEV_NB

        # Calculate the ratio 4jets/2jets in pp and pPb
        ratio_pp = int_sigma_4jets_pp/int_sigma_2jets_pp
        ratio_pPb = int_sigma_4jets_pPb/int_sigma_2jets_pPb

        plt.suptitle("Ratio of 4-jet to 2-jet production")
        ax.set_ylabel(r"Jet production ratio $R = \sigma_{\text{4-jets}}/\sigma_{\text{2-jets}}$")

        # Add Monte Carlo integration setup data
        ax.text(0.60, 0.04, rf"Monte Carlo: $N={N}$ with {n_bins} bins", color="gray", transform=ax.transAxes, fontsize=9.5, verticalalignment="bottom")
        
        plt.scatter(Delta_y, ratio_pp, color='k', label='pp', s=16)
        plt.scatter(Delta_y, ratio_pPb, color='red', label='pPb', s=16)
    
    # if DPS ratio in proton-lead and proton-proton collisions
    elif ratio_type == "dps":
        logger.info("Plotting DPS ratio in pp versus pPb")
        # Calculate the DPS ratio in nuclear vs pure proton collisions 
        ratio_dps = int_sigma_4jets_pPb/(A*int_sigma_4jets_pp)
        print(f"Mean DPS ratio: {np.mean(ratio_dps)}")

        plt.suptitle("4-jet DPS ratio of proton-lead versus proton-proton collisions")
        ax.set_ylabel(r"DPS ratio $R = 1/A \, \cdot \, \sigma_{pPb \rightarrow \text{4 jets}}/\sigma_{pp \rightarrow \text{4 jets}}$")

        # Add Monte Carlo integration setup data
        ax.text(0.60, 0.07, rf"Monte Carlo: $N={N}$ with {n_bins} bins", color="gray", transform=ax.transAxes, fontsize=9.5, verticalalignment="bottom")

        plt.scatter(Delta_y, ratio_dps, color='red', s=16)
        # Plot a horizontal line at R=1 to indicate no nuclear modifications
        ax.hlines(1, xmin=y_range_min, xmax=y_range_max, label="No nuclear modifications", linestyle='--', color='k')

    else: raise ValueError("Unknown ratio_type parameter...")

    plt.legend()
    plt.grid(alpha=0.4)
    plt.show()


def calculate_dps_ratio(pt_cut:float, N_total:int=10_000, N_injet:int=1_000_000) -> tuple:
    """
    Function calculating the (in-jet) DPS ratio in proton-lead and proton-proton collisions, based on the output of the Monte Carlo integration employed in double_dijet_sigma_total.
    This is an alternative approach to the above plotting in rapidity space, simply producing a float type number for comparison through an integration over the entire available phase space.
    Note that N_injet >> N_total, which is necessitated by the Heaviside condition. Ideally, the number of points satisfying the kinematical condition, used to evaluate the integral, should be approximately the number of N_total, to ensure similar Monte Carlo accuracy.  

    :param pt_cut: minimum transverse momentum of a jet pair, used as a cut in the cross section evaluation
    :param N_total: number of Monte Carlo samples to generate for each parameter in the integration scheme of the total cross section
    :param N_injet: number of Monte Carlo samples to generate for each parameter in the integration scheme of the in-jet cross section

    :returns: total and in-jet DPS ratios, normalised through the nucleon number of the nucleus A
    """

    # Initialise the empty ratio list, to average over as the final result
    dps_ratio_total = np.array([])
    dps_ratio_injet = np.array([])

    # Calculate the cross sections, total or injet-constrained
    logger.info(f"Running Monte Carlo integration for total DPS ratio calculations w/ N={N_total}, pt_cut={pt_cut} GeV...")
    sigma_pp_total, _ = double_dijet_sigma_total(N_total, "pp", pt_cut)
    sigma_pPb_total, _ = double_dijet_sigma_total(N_total, "pPb", pt_cut)
    dps_ratio_total = np.append(dps_ratio_total, sigma_pPb_total/(A*sigma_pp_total))

    logger.info(f"Running Monte Carlo integration for injet DPS ratio calculations w/ N={N_injet}, pt_cut={pt_cut} GeV...")
    sigma_pp_injet, _ = injet_double_dijet_sigma_total(N_injet, "pp", pt_cut, RADIUS)
    sigma_pPb_injet, _ = injet_double_dijet_sigma_total(N_injet, "pPb", pt_cut, RADIUS)
    dps_ratio_injet = np.append(dps_ratio_injet, sigma_pPb_injet/(A*sigma_pp_injet))
    
    return np.mean(dps_ratio_total, dtype=float), np.mean(dps_ratio_injet, dtype=float)


def plot_dps_ratio_vs_momentum_cut(pt_cut:list) -> None:
    """
    Plotting of the DPS ratio (sigma_pPb / (A*sigma_pp)) versus minimum transverse momentum cut, used to regulate the divergence in the cross section.

    :param pt_cut: list of momentum cuts to apply
    """

    # Initialise empty arrays
    ratio_ptcut_total = np.array([])
    ratio_ptcut_injet = np.array([])

    # Iterate over momentum cuts and calculate the DPS ratio for total and injet cross sections
    for cut in pt_cut:
        dps_ratio_total, dps_ratio_injet = calculate_dps_ratio(cut)
        ratio_ptcut_total = np.append(ratio_ptcut_total, dps_ratio_total)
        ratio_ptcut_injet = np.append(ratio_ptcut_injet, dps_ratio_injet)

    # Plotting
    fig, ax = plt.subplots(figsize=(9,6))
    plt.title("DPS ratio versus transverse momentum cut")
    ax.set_xlabel("Momentum cut [GeV]")
    ax.set_ylabel("DPS ratio")

    plt.scatter(pt_cut, ratio_ptcut_total, color='k', label='total')
    plt.scatter(pt_cut, ratio_ptcut_injet, color='red', label='in-jet')

    plt.legend()
    plt.grid(alpha=0.4)
    plt.show()


def plot_percentage_vs_momentum_cut(pt_cut:np.ndarray, N_total:int=10_000, N_injet:int=1_000_000, collision_type:str="pPb") -> None:
    """
    Plotting the relative percentage of events satisfying the injet kinematical condition, through the ratio in-jet/total, as a function of the applied minimal jet pair transverse momentum cut.

    :param pt_cut: minimum transverse momentum of a jet pair, used as a cut in the cross section evaluation
    :param N_total: number of Monte Carlo samples to generate for each parameter in the integration scheme of the total cross section
    :param N_injet: number of Monte Carlo samples to generate for each parameter in the integration scheme of the in-jet cross section
    :param collision_type: string representing the type of the interaction, i.e. either pp ("proton") or pA ("nuclear"); default setting corresponds to proton-proton collisions
    """ 
    
    # Initialise an empty array
    ratio =  np.array([])

    # Get cross sections for different pt cuts
    for cut in pt_cut:

        # Calculate the cross sections and the corresponding ratio
        logger.info(f"Running Monte Carlo integration for injet DPS ratio calculations w/ N={N_injet}, pt_cut={cut} GeV...")
        sigma_injet, _ = injet_double_dijet_sigma_total(N_injet, collision_type, cut, RADIUS)

        logger.info(f"Running Monte Carlo integration for total DPS ratio calculations w/ N={N_total}, pt_cut={cut} GeV...")
        sigma_total, _ = double_dijet_sigma_total(N_total, collision_type, cut)

        ratio = np.append(ratio, sigma_injet/sigma_total * 100) # in %
    
    # Plotting
    fig, ax = plt.subplots(figsize=(9,6))
    plt.title(f"Relative in-jet/total ratio versus transverse momentum cut in {collision_type} collisions")
    ax.set_xlabel("Momentum cut [GeV]")
    ax.set_ylabel("in-jet/total ratio [%]")

    plt.scatter(pt_cut, ratio, color='k')

    plt.grid(alpha=0.4)
    plt.show()


def plot_percentage_vs_radius(radii:np.ndarray, N_total:int=10_000, N_injet:int=1_000_000, collision_type:str="pPb") -> None:
    """
    Plotting the relative percentage of events satisfying the injet kinematical condition, through the ratio in-jet/total, as a function of the (parametrically small) radius of the jet cone.

    :param radii: array of cone radii
    :param N_total: number of Monte Carlo samples to generate for each parameter in the integration scheme of the total cross section
    :param N_injet: number of Monte Carlo samples to generate for each parameter in the integration scheme of the in-jet cross section
    :param collision_type: string representing the type of the interaction, i.e. either pp ("proton") or pA ("nuclear"); default setting corresponds to proton-proton collisions
    """

    # Define a parabola to fit to the data ~ hypothesis: the ratio grows with the square of the jet radius
    def square(R:np.ndarray, a:float, b:float) -> np.ndarray:
        return a*R**2 + b

    # Initialise an empty array
    ratio = np.array([])
    
    # Get cross sections for different radii
    for radius in radii:

        # Calculate the cross sections and the corresponding ratio @ PT_CUT momentum cut
        logger.info(f"Running Monte Carlo integration for injet DPS ratio calculations w/ N={N_injet}, R={radius} @ {PT_CUT} GeV...")
        sigma_injet, _ = injet_double_dijet_sigma_total(N_injet, collision_type, PT_CUT, radius)

        logger.info(f"Running Monte Carlo integration for total DPS ratio calculations w/ N={N_total} @ {PT_CUT} GeV...")
        sigma_total, _ = double_dijet_sigma_total(N_total, collision_type, PT_CUT)

        ratio = np.append(ratio, sigma_injet/sigma_total * 100) # in %

    # Fit the R^2 curve using SciPy
    popt, _ = curve_fit(square, radii, ratio)

    # Plotting
    fig, ax = plt.subplots(figsize=(9,6))
    plt.title(f"Relative in-jet/total ratio versus jet cone radius in {collision_type} collisions")
    ax.set_xlabel(r"Jet cone radius $R$")
    ax.set_ylabel("in-jet/total ratio [%]")

    plt.scatter(radii, ratio, color='k')
    plt.plot(radii, square(radii, *popt), color='r', label=f'fit a={popt[0]:.2f}, b={popt[1]:.2f}')

    plt.legend()
    plt.grid(alpha=0.4)
    plt.show()


def fetch_sigma(PT:float, N_injet:int=1_000_000, collision_type:str="pPb", radius:float=RADIUS, pt_cut:float=PT_CUT) -> float:
    """
    Function returning the differential injet cross section evluated with global parameters PT_CUT and RADIUS, defined for the purposes of multiprocessing.

    :param PT: array of transverse momenta as reconstructed by the experiment
    :param N_injet: number of Monte Carlo samples to generate for each parameter in the integration scheme of the in-jet cross section  
    :param collision_type: string representing the type of the interaction, i.e. either pp ("proton") or pA ("nuclear"); default setting corresponds to proton-proton collisions
    :param radius: value of the jet radius determined by experiment
     :param pt_cut: float value of the momentum cut used to regulate the divergences in the cross section

    :returns: differential injet cross section
    """
    sigma_dPT, _ = injet_double_dijet_sigma_dPT(PT, N_injet, collision_type, pt_cut, radius)
    return sigma_dPT

def fetch_sigma_radius(radius:float, PT:float=30.0, N_injet:int=1_000_000, collision_type:str="pPb", pt_cut:float=PT_CUT) -> float:
    """
    Function returning the differential injet cross section evluated with global parameters PT_CUT and RADIUS, defined for the purposes of multiprocessing.

    :param pt_cut: float value of the momentum cut used to regulate the divergences in the cross section
    :param PT: value of transverse jet momentum as reconstructed by the experiment (ensuring PT >> pt_cut)
    :param N_injet: number of Monte Carlo samples to generate for each parameter in the integration scheme of the in-jet cross section  
    :param collision_type: string representing the type of the interaction, i.e. either pp ("proton") or pA ("nuclear"); default setting corresponds to proton-proton collisions
    :param radius: value of the jet radius determined by experiment

    :returns: differential injet cross section
    """
    sigma_dPT, _ = injet_double_dijet_sigma_dPT(PT, N_injet, collision_type, pt_cut, radius)
    return sigma_dPT

def fetch_sigma_cut(pt_cut:float, radius:float=RADIUS, PT:float=30.0, N_injet:int=500_000, collision_type:str="pPb") -> float:
    """
    Function returning the differential injet cross section evluated with global parameters PT_CUT and RADIUS, defined for the purposes of multiprocessing.

    :param pt_cut: float value of the momentum cut used to regulate the divergences in the cross section
    :param PT: value of transverse jet momentum as reconstructed by the experiment (ensuring PT >> pt_cut)
    :param N_injet: number of Monte Carlo samples to generate for each parameter in the integration scheme of the in-jet cross section  
    :param collision_type: string representing the type of the interaction, i.e. either pp ("proton") or pA ("nuclear"); default setting corresponds to proton-proton collisions
    :param radius: value of the jet radius determined by experiment

    :returns: differential injet cross section
    """
    sigma_dPT, _ = injet_double_dijet_sigma_dPT(PT, N_injet, collision_type, pt_cut, radius)
    return sigma_dPT

def plot_differential_vs_momentum(PT:np.ndarray) -> None:
    """
    Plotting the differential injet cross section vs the total transverse momentum P_T of the reconstructed jet, optimised to use multiple CPUs.

    :param PT: array of transverse momenta as reconstructed by the experiment
    """
    
    # Define a model inverse quartic function for fitting
    def quartic(x:np.ndarray, a:float, b:float) -> np.ndarray: 
        return a/(x**4) + b
    
    def analytical(x:np.ndarray, a:float, b:float) -> np.ndarray:
        return a*(1/(PT_CUT**2 * x**4) + 6/(PT_CUT * x**5) + 12/x**6 * np.log(x/PT_CUT) - 7/x**6) + b
    
    def square(x:np.ndarray, a:float, b:float) -> np.ndarray:
        return a/(x**2) + b

    # Split the task across available CPUs; some optimisation needed to manage WSL's memory cap
    with Pool(processes=(cpu_count()-1), maxtasksperchild=1) as pool:
        sigma_arr = np.array(pool.map(fetch_sigma, PT))

    # Try fitting with an inverse quartic function and expected analyical behaviour with scipy
    popt, _ = curve_fit(quartic, PT, sigma_arr)
    popt_analytical, _ = curve_fit(analytical, PT, sigma_arr)
    popt_square, _ = curve_fit(square, PT, sigma_arr)

    # Plotting
    fig, ax = plt.subplots(figsize=(9,6))
    plt.title(f"Differential in-jet DPS cross section vs reconstructed jet momentum @ pt_cut={PT_CUT} GeV")
    ax.set_xlabel(r"Jet transverse momentum $P_T$ [GeV]")
    ax.set_ylabel(r"$1/P_T$ $\cdot$ $d\sigma_{\text{in-jet}}/dP_T$ [GeV$^{-4}$]")
    
    plt.scatter(PT, sigma_arr, color='k', label="Monte Carlo")
    plt.plot(PT, quartic(PT, *popt), color='red', label=f"fit quartic a={popt[0]:.2}, b={popt[1]:.2}")
    #plt.plot(PT, analytical(PT, *popt_analytical), color='green', label=f"fit analytical a={popt_analytical[0]:.2}, b={popt_analytical[1]:.2}")
    plt.plot(PT, square(PT, *popt_square), color='blue', label=f"fit square a={popt_square[0]:.2}, b={popt[1]:.2}")

    plt.legend()
    plt.grid(alpha=0.4)
    plt.show()


def plot_dps_vs_sps(PT:np.ndarray) -> None:
    """
    Function plotting the ratio between the double parton scattering production of 4 jets, with the applied kinematical constraint, and the sinlge parton scattering production of two jets, in various jet transverse momentum intervals.

    :param PT: array of transverse momenta as reconstructed by the experiment
    """
    
    # Get the DPS and SPS differential cross sections
    # Split the task across available CPUs; some optimisation needed to manage WSL's memory cap
    with Pool(processes=(cpu_count()-1), maxtasksperchild=1) as pool:
        dps_sigma_arr = np.array(pool.map(fetch_sigma, PT))
        sps_sigma_arr = np.array(pool.map(dijet_sigma_dpt, PT))

    # Compute the ratio
    ratio = dps_sigma_arr / sps_sigma_arr

    # Plotting
    fig, ax = plt.subplots(figsize=(9,6))
    plt.title(f"Ratio of differential DPS/SPS cross sections vs reconstructed jet momentum @ pt_cut={PT_CUT}")
    ax.set_xlabel(r"Jet transverse momentum $P_T$ [GeV]")
    ax.set_ylabel(r"Sigma ratio $R(p_T)$")
    
    plt.scatter(PT, ratio, color='k', label="Monte Carlo")

    plt.legend()
    plt.grid(alpha=0.4)
    plt.show() 


def plot_differential_vs_radius(radius:np.ndarray) -> None:
    """
    Plotting the 1/PT dsigma/dPT differential injet cross section versus jet radius.

    :param: radius: array of radii for which to plot the cross section
    """
    
    # Define expected behaviour function for fitting; a defines the prefactors stemming from the PDF-dependence
    def func(r:np.ndarray, a:float, b:float) -> np.ndarray: 
        return a*r + b*r**2

    # Use multiprocessing to evaluate the cross sections for several values of the momentum cut
    with Pool(processes=(cpu_count() - 1), maxtasksperchild=1) as pool:
        dPT_sigma_arr = np.array(pool.map(fetch_sigma_radius, radius))

    # Try fitting with an inverse square function with scipy
    popt_square, _ = curve_fit(func, radius, dPT_sigma_arr)

    # Plotting
    fig, ax = plt.subplots(figsize=(9,6))
    plt.title(f"Differential in-jet DPS cross section versus jet radius @ PT=30 GeV")
    ax.set_xlabel(r"Radius $R$")
    ax.set_ylabel(r"$1/P_T$ $\cdot$ $d\sigma_{\text{in-jet}}/dP_T$ [GeV$^{-4}$]")
    
    plt.scatter(radius, dPT_sigma_arr, color='k', label="Monte Carlo")
    plt.plot(radius, func(radius, *popt_square), color='blue', label=f"fit a={popt_square[0]:.2}, b={popt_square[1]:.2}")

    print(popt_square[0]/((-1)*4*popt_square[1]))

    plt.legend()
    plt.grid(alpha=0.4)
    plt.show() 


def plot_differential_vs_cut(pt_cut:np.ndarray) -> None:
    """
    Plotting the differential injet cross section versus the momentum cut at PT defined by the fetching function.

    :param pt_cut: array of momentum cuts
    """

    PT = 30.0

    def square(x:np.ndarray, a:float, b:float):
        return a/(PT**4*x**2) + b

    def func(x:np.ndarray, a:float, b:float):
        return a*(1/(PT**4 * x**2) + 6/(PT**5 * x) + 12/PT**6 * np.log(PT/x) - 7/PT**6) + b

    with Pool(processes=(cpu_count() - 1), maxtasksperchild=1) as pool:
        dPT_sigma_arr = np.array(pool.map(fetch_sigma_cut, pt_cut))

    # Try fitting
    popt_func, _ = curve_fit(func, pt_cut, dPT_sigma_arr)
    popt_square, _ = curve_fit(square, pt_cut, dPT_sigma_arr)

    # Plotting
    fig, ax = plt.subplots(figsize=(9,6))
    plt.title(f"Differential in-jet DPS cross section versus transverse jet momentum cut @ PT=30 GeV")
    ax.set_xlabel(r"Transverse momentum cut $p_{T_{\text{cut}}}$")
    ax.set_ylabel(r"$1/P_T$ $\cdot$ $d\sigma_{\text{in-jet}}/dP_T$ [GeV$^{-4}$]")
    
    plt.scatter(pt_cut, dPT_sigma_arr, color='k', label="Monte Carlo")
    plt.plot(pt_cut, func(pt_cut, *popt_func), color='green', label=f"fit analytical a={popt_func[0]:.2}, b={popt_func[1]:.2}")
    plt.plot(pt_cut, square(pt_cut, *popt_square), color='red', label=f"fit square a={popt_square[0]:.2}, b={popt_square[1]:.2}")

    plt.legend()
    plt.grid(alpha=0.4)
    plt.show() 




# --- Main execution ---
def main() -> None:
    pass


if __name__ == "__main__":
    main()
# ---------------------