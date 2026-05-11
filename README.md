## Description
The following is an implementation of the calculation of the contribution of double parton scattering (DPS) to the leading-order cross section for four-jet production in proton-proton and proton-nucleus collisions. The main flow of the programme is to first compute the single parton scattering cross section for production of two jets, to be treated as two separate hard processes, in the implementation of the "pocket formula" for DPS into 4 jets. Moreover, in order to estimate the extent of this effect in proton-nucleus collisions, the functions used in the proton-proton case are characterised by a degree of generality, with a possible input representing the nucleon interacting with the proton in the course of the collision.

The mathematical treatment employed in this package follows from ~ J. F. Owens' "Large momentum-transfer production of direct photons, jets, and particles", where the cross section for the process $AB \rightarrow \text{2 jets}$ takes the following form

$$ \frac{d\sigma_{AB \rightarrow \text{jet}_1 + \text{jet}_2}}{dy_1dy_2dp_T^2} = \sum_{ab} x_a f_{a/A}(x_a) \, x_b f_{b/B}(x_b) \frac{d\sigma}{d\hat{t}}(ab \rightarrow 1 2),$$

where $f_{a,b/A,B}$ represent the distribution functions for partons $a$ or $b$ within hadrons $A$ or $B$, respectively, and $d\sigma/d\hat{t}$ is the partonic differential cross section for $2\rightarrow 2$ scattering, at leading order in QCD. As such, for the collision of two protons, we choose $A=B=p$, however, if one considers nuclear collisions, depending on the model, collisions including neutrons also need to be accounted for.

The particular modelling chosen for the proton-nucleus collision follows ~ H-S. Shao's "Probing impact-parameter dependent nuclear parton densities from double parton scatterings in heavy-ion collisions". As a first approach, no spatial dependence of the nuclear parton distribution functions is assumed. Moreover, a hard-sphere framework is chosen for the modelling of the nucleon density. 

$$ \frac{d\sigma^{\text{DPS}}_{Ap\rightarrow 4 \text{ jets}}}{dy_1dy_2dy_3dy_4dp_{T_1}^2dp_{T_2}^2} =  \sum_{N_1^A, N_2^A} \frac{d\sigma_{N_1^A p\rightarrow 2 \text{ jets}}}{dy_1dy_2dp_{T_1}^2} \frac{d\sigma_{N_2^Ap\rightarrow 2 \text{ jets}}}{dy_3dy_4dp_{T_2}^2} \left(\frac{\delta_{N_1^A N_2^A}}{\sigma_{\text{eff}}} + \frac{A-1}{A}\frac{9}{8\pi R_A^2}\right)$$

The sum in the above expression are performed over all nucleons in the nucleus $A$. The term for which $N_1^A = N_2^A$ represent the 4-parton interaction in which two partons originate from a single nucleon inside the nucleus, whereas the one where $N_1^A \neq N_2^A$ signifies the situation where they stem from a combination of two individual nucleons, either proton-proton, neutron-neutron, proton-neutron, or neutron-proton.

The specific proton-nucleus interaction treated in this implementation is the proton-lead (A=208) collision, utilising CT18ANLO free proton PDFs and EPPS21 nuclear PDFs, encompassed in the pre-compiled LHAPDF Python library.

## Package structure
- `process_vars.py`: List of global variables, defining either the collision process, or the characteristics of the jets. Some of these variables, such as the min/max jet transverse momentum or max rapidity, are crucial for setting proper integration limits in the cross section computations. 

- `load_pdf.py`: Module loading the parton distribution functions using the LHAPDF library and the specific set determined in `process_vars`.

- `alpha_s.py`: Encoding the one-loop running of the strong coupling constant. The default convention is to use 4 quark flavours and a Lambda4 scale parameter, following the convention set by the NNPDF4.0 PDF set @ NNLO. However, if the PDF set specified in `process_vars` includes a fitting of the coupling running, the function will use the values instead, in order to maintain consistency within the global fit.

- `partonic_sigma.py`: Calculating the $\hat{s}^2/\pi \alpha_s^2 \cdot d\sigma/d\hat{t}$ of parton-level 2-to-2 scatterings at leading order in  perturbative QCD, relevant for the formation of the dijets. The reason for such a parametrisation of the output quantities is to simplify the calculations by including the pre-factor once at a later stage, instead of at the level of each partonic cross section.

- `dijet_sigma.py`: Calculating the single parton scattering contribution to the differential cross section for $pA \rightarrow 2\text{ jets}$, i.e. $d\sigma/dy_1 dy_2 dp_T^2$. Included is also an appropriate integration to express the differential cross section in terms of rapidity differences, for a range of transverse jet momentum values.

- `double_dijet_sigma.py`: Using the factorised ansatz, the double parton scattering contribution is calculated, using pre-computed results from the single parton scattering in `dijet_sigma`, resulting in the quantity $d\sigma/dy_1dy_2dy_3dy_4dp_{T_1}^2dp_{T_2}^2$ Moreover, the result is integrated using a Monte Carlo scheme, to express the differential cross section in terms of maximal rapidity separation between the two most remote jets in the configuration, i.e. $d\sigma/d\Delta y$.

- `jet_overlap_sigma.py`: Building on the differential DPS cross section, the function `injet_double_sigma_total` calculates the total cross section corresponding to finding jet 3 within jet 1. This is done by enforcing a Heaviside constraint, defining the original jet cone of radius $R$, and restoring angular-dependence of the expression by averaging over the azimuthal angles of the two involved jets. Adapting the previous approaches, the full integration is performed using a Monte Carlo method.

- `main.py`: Main execution loop of the programme: included are functions visualising the results of the cross section calculation, based on various criteria. `plot_sigma()` produces a log-linear plot of the differential cross section, expressed in terms of rapidity separation for proton-proton collisions, in both DPS and SPS interactions. On the other hand, `plot_ratio()` calculates the ratio of 4-jet to 2-jet production in proton-proton and proton-lead collisions, plotting it on a common scale, quantifying the appearance of DPS effects in a nuclear environment. It also evaluates the ratio of 4-jet production in proton-lead versus proton-proton interactions, establishing the extent of nuclear modifications in the same process.