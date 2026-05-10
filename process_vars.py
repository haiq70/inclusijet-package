import lhapdf

SQRT_S = 5000 # Centre-of-mass energy in GeV
S = SQRT_S**2 # in GeV^2
SIGMA_EFF = 11.3 * 2.56819 # Effective pp cross-section in GeV^{-2} (15 mb)

# Convert GeV^{-2} to nb
CONV_GEV_NB = 3.89379e5

# Global Monte Carlo numbering scheme to refer to the partons
MC_IDX = {
    "bbar": -5, "cbar": -4, "sbar": -3, "ubar": -2, "dbar": -1,
    "d": 1, "u": 2, "s": 3, "c": 4, "b": 5, 
    "g": 21
}

# Set up dictionary mapping for the action of isospin symmetry: takes desired neutron pdf parton input and maps to the corresponding proton pdf parton
ISOSPIN_MAP = {"u": "d", "d": "u", "ubar": "dbar", "dbar": "ubar", "s": "s", "sbar": "sbar", "g": "g"}

# Relevant partons to take into account in the calculation, evaluating PDFs/FFs
QUARKS = ["u", "d", "s"]
ANTIQUARKS = ["ubar", "dbar", "sbar"] 
FLAVOURS = QUARKS + ANTIQUARKS + ["g"]


# --- Jets ---
# Assume a range of transverse momenta
PT_MIN = 20
PT_MAX = 40

# Jet radius in angular space
RADIUS = 0.4

# Assume same rapidity range for the jets |y| < Y_MAX
Y_MAX = 4.7
# ---------

# --- Nucleus ---
# Nucleon and atomic number for lead
A = 208
Z = 82

# Nuclear radius
# 1 fm = 5.07 GeV^{-1}
RA =  6.624 * 5.07 # GeV^{-1}
# ---------------

# Specification of PDF sets
PDF_PROTON_TITLE = "CT18ANLO"
PDF_PROTON = lhapdf.mkPDF(PDF_PROTON_TITLE, 0)

PDF_NUCLEUS_TITLE = "EPPS21nlo_CT18Anlo_Pb208"
PDF_NUCLEUS = lhapdf.mkPDF(PDF_NUCLEUS_TITLE, 0)