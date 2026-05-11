# --- Partonic cross sections ---
# Quarks denoted as q, antiquarks as b, gluons as g

def sigqiqjqiqj(shat:float, that:float, uhat:float) -> float:
    """
    Rescaled parton-level differential cross section for the qq'->qq' channel.

    :param shat: partonic Mandelstam s variable.
    :param that: partonic Mandelstam t variable.
    :param uhat: partonic Mandelstam u variable.

    :returns: rescaled differential cross section
    """
    return 4/9 * (shat**2 + uhat**2)/that**2

def sigqiqiqiqi(shat:float, that:float, uhat:float) -> float:
    """
    Rescaled parton-level differential cross section for the qq->qq channel.

    :param shat: partonic Mandelstam s variable.
    :param that: partonic Mandelstam t variable.
    :param uhat: partonic Mandelstam u variable.

    :returns: rescaled differential cross section
    """
    return 4/9 * ((shat**2 + uhat**2)/that**2 + (shat**2 + that**2)/uhat**2) - 8/27 * shat**2/(that * uhat)

def sigqibiqjbj(shat:float, that:float, uhat:float) -> float:
    """
    Rescaled parton-level differential cross section for the q qbar->q' qbar' channel.

    :param shat: partonic Mandelstam s variable.
    :param that: partonic Mandelstam t variable.
    :param uhat: partonic Mandelstam u variable.

    :returns: rescaled differential cross section
    """
    return 4/9 * (that**2 + uhat**2)/shat**2


def sigqibiqibi(shat:float, that:float, uhat:float) -> float:
    """
    Rescaled parton-level differential cross section for the q qbar->q qbar channel.

    :param shat: partonic Mandelstam s variable.
    :param that: partonic Mandelstam t variable.
    :param uhat: partonic Mandelstam u variable.

    :returns: rescaled differential cross section
    """
    return 4/9 * ((shat**2 + uhat**2)/that**2 + (uhat**2 + that**2)/shat**2) - 8/27 * uhat**2/(shat*that)

def siggqigqi(shat:float, that:float, uhat:float) -> float:
    """
    Rescaled parton-level differential cross section for the gq->gq channel.

    :param shat: partonic Mandelstam s variable.
    :param that: partonic Mandelstam t variable.
    :param uhat: partonic Mandelstam u variable.

    :returns: rescaled differential cross section
    """
    return -4/9 * (shat/uhat + uhat/shat) + (shat**2 + uhat**2)/that**2

def sigqibigg(shat:float, that:float, uhat:float) -> float:
    """
    Rescaled parton-level differential cross section for the q qbar->gg channel.

    :param shat: partonic Mandelstam s variable.
    :param that: partonic Mandelstam t variable.
    :param uhat: partonic Mandelstam u variable.

    :returns: rescaled differential cross section
    """
    return 32/27 * (that/uhat + uhat/that) - 8/3 * (that**2+uhat**2)/shat**2

def sigggqibi(shat:float, that:float, uhat:float) -> float:
    """
    Rescaled parton-level differential cross section for the gg->q qbar channel.

    :param shat: partonic Mandelstam s variable.
    :param that: partonic Mandelstam t variable.
    :param uhat: partonic Mandelstam u variable.

    :returns: rescaled differential cross section
    """
    return 1/6 * (that/uhat + uhat/that) - 3/8 * (that**2 + uhat**2)/shat**2

def siggggg(shat:float, that:float, uhat:float) -> float:
    """
    Rescaled parton-level differential cross section for the gg->gg channel.

    :param shat: partonic Mandelstam s variable.
    :param that: partonic Mandelstam t variable.
    :param uhat: partonic Mandelstam u variable.

    :returns: rescaled differential cross section
    """
    return 9/2 * (3 - that*uhat/shat**2 - shat*uhat/that**2 - shat*that/uhat**2)

# --------------------------------------------