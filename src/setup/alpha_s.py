import numpy as np

from setup.process_vars import PDF_PROTON_TITLE, PDF_PROTON, PDF_NUCLEUS_TITLE

# --- Extract running alpha_s coupling value ---
def alpha_s(Q2:float, nf:int=4, Lambda:float=0.342207) -> float:
    """
    Implementing the one-loop running of the strong coupling as a function of Q^2.
    
    :param Q2: characteristic squared energy scale
    :param nf: number of quark flavours relevant at a given Q^2 (default = 4)
    :param Lambda: flavour-dependent energy scale where alpha_s becomes strong as Q^2 is decreased (default nf=4: 0.342207 GeV in NNPDF4.0 @ NNLO)
    
    :returns: alpha_s value for a given Q^2
    """

    Q = np.sqrt(Q2)

    # Handle negative log values
    if Q2 <= Lambda**2: return 0.0

    # If the PDF set contains alpha_s evolution (e.g. CT18), use it instead of one-loop RG evolution
    if PDF_PROTON_TITLE == "CT18ANLO" or PDF_NUCLEUS_TITLE == "EPPS21nlo_CT18Anlo_Pb208":
        # Note that the alpha_s evolution is the identical for the nuclear and proton baseline sets
        return PDF_PROTON.alphasQ(Q)

    # Note that the default Lambda is to be changed, depending on the used PDF set and the number of effective flavours
    return 12*np.pi/(33 - 2*nf)/np.log(Q2/Lambda**2)
# ---------------------------------------