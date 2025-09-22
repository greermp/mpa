from models import *
from simulation.doe import create_doe
import math

import numpy as np
debug = True

def get_sensor_footprint(alt_ft, sensor):
    """
    Compute sensor footprint based on altitude and FOV.

    Parameters
    ----------
    alt_ft : float
        Altitude in ft.
    sensor : EO_Sensor
        Sensor with fov_deg_w and fov_deg_l

    Returns
    -------
    float
        Cost in millions.
    """
    footprint_length = alt_ft * 0.0003048 *2 * math.tan(sensor.fov_deg_w/2*math.pi/180)
    footprint_width  = alt_ft * 0.0003048 *2 * math.tan(sensor.fov_deg_l/2*math.pi/180)
    return footprint_length * footprint_width 

def mach_to_mps(mach, alt_ft):
    """
    Compute meters/second based on altitude and mach
    
    a = 340.3 m/s * sqrt(T/288.15)
    
    ISA temp lapse: 15C at sea level, -2 degreesC per 1000ft )
    Kelvin: T(h) = 288.15 - 0.0065 x h(meters) 

    Parameters
    ----------
    mach : float
        Mach number
    alt_ft : float
        Altitude in ft.

    Returns
    -------
    float
        Speed in m/s.
    """
    alt_m    = alt_ft * 0.3048
    t_kelvin = 288.15 - 0.0065 * alt_m
    a        = 340.3 * math.sqrt(t_kelvin/288.15) # a: Speed of Sound (m/s)
    v_ms     = mach * a
    return v_ms
    
def time_to_fly_xkm(mps, xkm):
    distance_m = xkm * 1000  
    # time in seconds
    time_sec = distance_m / mps  
    # convert to hours
    return time_sec / 3600  
    
def calc_area_sanitized(mps, endurance, footprint, ingress_dist_km, x=100, y=100):
    # mps = mach_to_mps(mach, alt_ft)
    if (debug) : print("Starting endurance: ", endurance )
    ingress_time = time_to_fly_xkm(mps, ingress_dist_km)
    endurance -= ingress_time
    if (debug) : print("After ingress: ", endurance )
    
    x_patrol_time = time_to_fly_xkm(mps, x)
    if (debug) : print("Patrol time for 1 strip (x): ", x_patrol_time )
    y_complete = 0
    while y_complete < 100:
        if endurance >=  x_patrol_time:
            endurance -= x_patrol_time
            # y_complete = 
    exit

def main():
    ingress_dist_km = 100
    ### Factors
    # Speed (Mach) from 0.4 to 0.9, e.g., 6 levels
    machs = np.linspace(0.4, 0.9, 6)
    # Altitude from 5k to 25k feet
    altitudes_ft = np.linspace(5000, 25000, 5)
    
    # Assume a spacing between sweeps of footprint * (1-sweep_overlap) 
    sweep_overlap = 0.2

    doe = create_doe(machs, altitudes_ft)
    
    # compute velocity, cost, endurance, sensor footprint
    doe['speed_mps'] = doe.apply(
        lambda row: mach_to_mps(row.speed_mach, row.alt_ft), axis=1)
    
    # Keep cost function in kft...
    doe['cost'] = doe.apply(
        lambda row: cost_model(row.speed_mach, row.alt_ft/1000, row.sensor), axis=1)
    # Keep endurance function in kft...
    doe['endurance_hr'] = doe.apply(
        lambda row: endurance_model(row.speed_mach, row.alt_ft/1000), axis=1)
    
    doe['footprint_km2'] = doe.apply(
        lambda row: get_sensor_footprint(row.alt_ft, row.sensor), axis=1)

    if (False):
        for idx, row in doe.iterrows():
            # sanitized = calc_area_sanitized(row["speed_mps"], row["endurance_hr"], row["footprint_km2"], ingress_dist_km)
            sanitized = calc_area_sanitized(row["speed_mps"], row["endurance_hr"], row["footprint_km2"], ingress_dist_km)
            doe.loc[idx, "area_sanitized_km2"] = sanitized
            
            print(doe.head(1))
            break




if __name__ == "__main__":
    main()