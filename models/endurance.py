def endurance_model(mach, alt_kft):
    """
    Compute aircraft endurance in hours given Mach number and altitude.

    Parameters
    ----------
    mach : float
        Mach number
    altitude : float
        Altitude in kft.

    Returns
    -------
    float
        Endurance in hours.
    """
    return -18.75 * mach**2 + 8.0893 * mach + 0.01 * alt_kft**2 + 0.05 * alt_kft + 9.2105