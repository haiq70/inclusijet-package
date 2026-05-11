import numpy as np
import matplotlib.pyplot as plt
import logging

from dijet_sigma import dijet_sigma_dDelta_y
from double_dijet_sigma import double_dijet_sigma_Delta_y_max, double_dijet_sigma_total
from process_vars import Y_MAX, PT_MIN, PT_MAX, SQRT_S, PDF_PROTON_TITLE, PDF_NUCLEUS_TITLE, CONV_GEV_NB, A

plt.rcParams.update({
    "font.family":          "serif",
    "font.serif":           ["DejaVu Serif"],
    "xtick.minor.visible":  True,
    "ytick.minor.visible":  True,

})

# Set up logging info for programme execution
logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


def plot_sigma(Delta_y:np.ndarray, N:int=60000, n_bins:int=12, y_range_min:float=4.0, y_range_max:float=9.3, type:str="pp") -> None:
    """
    Visualising the results of the cross section for pp -> 2 jets and DPS pp -> 4 jets for a range of rapidity differences.

    :param Delta_y: rapidity differences between two jets
    :param N: Monte Carlo samples to generate
    :param n_bins: number of histogram bins
    :param y_range_min: minimum of the range of rapidity distances to plot over (x-values)
    :param y_range_max: maximum of the range of rapidity distances to plot over (x-values)
    :param type: string representing the type of the interaction, i.e. either pp ("proton") or pA ("nuclear"); default setting corresponds to proton-proton collisions
    """ 

    logger.info("Plotting SPS vs DPS sigma in pp...")

    # Run SPS
    def sps() -> np.ndarray:
        # Initialise a list of results for dsigma/dDelta_y
        dsigma = np.array([dijet_sigma_dDelta_y(i, type) for i in Delta_y]) * CONV_GEV_NB
        return dsigma

    # Run DPS
    def dps() -> tuple:
        logger.info(f"Running Monte Carlo integration for DPS w/ N={N} and n_bins = {n_bins}")
        res, bin_centres = double_dijet_sigma_Delta_y_max(N, n_bins, y_range_min, y_range_max, type)
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


def calculate_dps_ratio(N:int=30000):
    """
    Function calculating the DPS ratio in proton-lead and proton-proton collisions, based on the output of the Monte Carlo integration employed in double_dijet_sigma_total.
    This is an alternative approach to the above plotting in rapidity space, simply producing a float type number for comparison through an integration over the entire available phase space.

    :param N: number of Monte Carlo samples to generate for each parameter in the integration scheme

    :returns: DPS ratio, normalised through the nucleon number of the nucleus A
    """
    logger.info(f"Running Monte Carlo integration for DPS ratio calculations w/ N={N}...")

    sigma_pp, _ = double_dijet_sigma_total(N, type="pp")
    sigma_pPb, _ = double_dijet_sigma_total(N, type="pPb")

    return sigma_pPb/(A*sigma_pp)


# --- Main execution ---
def main() -> None:
    
    # y_range_min = 4.0
    # y_range_max = 9.3
    # type = "pp"

    # # Range of rapidity differences
    # Delta_y = np.arange(y_range_min, y_range_max, 0.1, dtype=float)

    # Plot the cross section of SPS and DPS in pp collisions
    # plot_sigma(Delta_y, type=type)

    # Plot jets or dps ratio
    #plot_ratio(ratio_type="dps")

    # Calculate the DPS ratio
    dps_ratio = calculate_dps_ratio()
    print(dps_ratio)


if __name__ == "__main__":
    main()
# ---------------------