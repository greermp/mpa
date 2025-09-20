# models.py

def endurance(mach, alt_kft):
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

# Sensor Cost in Millions
get_sensor_cost = {
    "EO/IR Sensor 1" : 50000/1e6,
    "EO/IR Sensor 2" : 1e6/1e6,
    "EO/IR Sensor 3" : 10e6/1e6
}

def cost(mach, alt_kft, sensor):
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
    try:
        sensor_cost = get_sensor_cost[sensor]
    except KeyError:
        raise ValueError(f"Unknown sensor: {sensor}")
    return 50 * mach**2 - 35 * mach + 0.03 * alt_kft**2 - 0.02 * alt_kft**2 + 11 + sensor_cost
