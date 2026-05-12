from setup.process_vars import MC_IDX, PDF_PROTON, PDF_NUCLEUS

# --- Fast PDF loading ---
def pdf_proton(parton:str, x:float, Q:float):
    """
    Return the proton PDF for a given parton at specified x, Q^2.

    :param parton: name of the parton, consistent with the Monte Carlo scheme
    :param x: momentum fraction variable
    :param Q: characteristic energy scale

    :returns: proton PDF
    """
    return PDF_PROTON.xfxQ(MC_IDX[parton], x, Q) / x


def pdf_nucleus(parton:str, x:float, Q:float):
    """
    Return the nuclear PDF for a given parton at specified x, Q^2.

    :param parton: name of the parton, consistent with the Monte Carlo scheme
    :param x: momentum fraction variable
    :param Q: characteristic energy scale

    :returns: nuclear PDF
    """
    return PDF_NUCLEUS.xfxQ(MC_IDX[parton], x, Q) / x
# ----------------------------