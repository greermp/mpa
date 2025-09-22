from models import *
from simulation.doe import create_doe
import math

import numpy as np
debug = False

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
    
def hours_to_fly_xkm(speed_mps, xkm):
    distance_m = xkm * 1000  
    # time in seconds
    time_sec = distance_m / speed_mps  
    # convert to hours
    return time_sec / 3600  
    
def calc_area_sanitized(speed_mps, endurance, footprint, ingress_dist_km, fuel_reserve_hrs, aoi_width=100, aoi_length=100):
    if debug : print(f"Starting endurance: {endurance}" )
    
    patrol_area = aoi_width*aoi_length
    ingress_time = hours_to_fly_xkm(speed_mps, ingress_dist_km)
    endurance -= ingress_time
    if debug : print(f"After ingress: {endurance}")
    endurance -= fuel_reserve_hrs
    if debug : print(f"After fuel reserve: {endurance}")
    
    x_patrol_time = hours_to_fly_xkm(speed_mps, aoi_width)
    if debug : print(f"Time for 1 strip (aoi_width): {x_patrol_time}" )
    area_sanitized = 0
    swath_width = math.sqrt(footprint) # FOV_Width
    
    distance_from_base_y = ingress_dist_km + (swath_width/2) # Start in center of swath
    
    while area_sanitized < patrol_area:
        rtb_time = hours_to_fly_xkm(speed_mps, distance_from_base_y)
        # if debug : print(f"{distance_from_base_y} km from base ({rtb_time}) hrs" )
        
        #TODO: Account for moving east/west, and changing distance from base
        if endurance <  x_patrol_time+rtb_time: # If we have enough TOS to complete this strip
            print("\nBingo fuel, RTB")
            if debug: print(f"Endurance remaining: {endurance}. Distance from base {distance_from_base_y}km.  Time for 1 more patrol and RTB:{hours_to_fly_xkm(speed_mps,aoi_width+distance_from_base_y)} ")
            return area_sanitized
        
        
        endurance                -= x_patrol_time            # deduct "fuel"
        area_sanitized           += swath_width*aoi_width    # record area sanitized
        distance_from_base_y     += swath_width              # move aircraft away from base
        progress_percent          = (area_sanitized / patrol_area) * 100
        print(f"\rPatrol progress: {progress_percent:.1f}% complete", end="", flush=True)
    return patrol_area
            
            
def main():
    ingress_dist_km  = 100
    fuel_reserve_hrs =   0.5
    ### Factors
    # Speed (Mach) from 0.4 to 0.9, e.g., 6 levels
    machs = np.linspace(0.4, 0.9, 6)
    # Altitude from 5k to 25k feet
    altitudes_ft = np.linspace(5000, 25000, 5)

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

    for idx, row in doe.iterrows():
        sanitized = calc_area_sanitized(row["speed_mps"],
                                        row["endurance_hr"],
                                        row["footprint_km2"],
                                        ingress_dist_km, fuel_reserve_hrs)
        doe.loc[idx, "area_sanitized_km2"] = sanitized
    
    doe.to_csv("mpa.csv", columns= [col for col in doe.columns if col != "sensor"] ,index=False)





if __name__ == "__main__":
    main()