# Sensor Cost in Millions

def cost_model(mach, alt_kft, sensor):
    """
    Compute aircraft cost (MM) given Mach number and altitude.

    Parameters
    ----------
    mach : float
        Mach number
    altitude : float
        Altitude in kft.

    Returns
    -------
    float
        Cost in millions.
    """
    # try:
    #     sensor_cost = sensor.co
    # except KeyError:
    #     raise ValueError(f"Unknown sensor: {sensor}")
    
    return 50 * mach**2 - 35 * mach + 0.03 * alt_kft**2 - 0.02 * alt_kft**2 + 11 + sensor.cost
