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

def percent_patrol_complete(gridsize_km2, area_sanitized):
    return area_sanitized/gridsize_km2*100

def distance_to_base(aircraft_x, aircraft_y, ingress_dist):
    """
    Compute distance to base (-ingres_dist, 0) to aircraft position
              AOI
             _ _ _
    (0,0)    _ _ _
    origin-> _ _ _ <--end first patrol leg (3,0)
            
            
            
            * Base (0, -ingress_dist)

    Parameters
    ----------
    aircraft_x (float): x-position within AOI (km)
    aircraft_y (float): y-position within AOI (km)
    ingress_dist (float): ingress distance from base to AOI entry (km)

    Returns
    -------
    float
        distance from base (km)
    """
    
    
def calc_area_sanitized(
    speed_mps,
    endurance,
    footprint,
    ingress_dist_km,
    fuel_reserve_hrs=0,
    aoi_width=100,
    aoi_length=100) -> float:
    """
    Simulate an MPA patrol over an AOI and compute the sanitized area before returning to base.

    The aircraft starts at the southern edge of the AOI (0,0) after ingress and flies
    strips northward, deducting endurance for each strip and the return trip to base.

    Base coordinates are at (0, -ingress_dist_km), directly south of AOI origin.

    AOI dimensions:
        Width  = aoi_width  (km, east-west)
        Length = aoi_length (km, north-south)

    Sensor footprint is assumed square; swath width = sqrt(footprint).

    Coordinates schematic (for ingress_dist_km = 100 km, aoi_width = 100 km):

          y (north)
          ↑
          |
          |       * Aircraft position (aircraft_x, aircraft_y)
          |       * After first strip, aircraft_y ≈ 0
          |
          * AOI origin (0,0) - southern edge of AOI
          |
          |
          * Base at (0, -ingress_dist_km)
          |
          +----------------------------→ x (east)

    Args:
        speed_mps (float): True airspeed of the aircraft in meters/second.
        endurance (float): Total available flight time (hours) including reserves.
        footprint (float): Sensor footprint area in km².
        ingress_dist_km (float): Distance from base to AOI southern edge (km).
        fuel_reserve_hrs (float): Minimum reserve endurance to retain (hours).
        aoi_width (float, optional): Width of AOI (km, default 100 km).
        aoi_length (float, optional): Length of AOI (km, default 100 km).

    Returns:
        float: Total area sanitized (km²) before hitting Bingo fuel and returning to base.
    """
    if debug : print(f"Starting endurance: {endurance}" )
    aircraft_x       = 0 # 0 representes W side of the grid, 1 the E side
    area_sanitized = 0
    patrol_area    = aoi_width*aoi_length
    x_patrol_time  = hours_to_fly_xkm(speed_mps, aoi_width) # time to sanitize from position (0, n) to (aoi_width, n)
    endurance     -= fuel_reserve_hrs # Subtract reserve fuel from `endurance`
    
    endurance -= hours_to_fly_xkm(speed_mps, ingress_dist_km) # Subtract time to ingress from `endurance`
    sensor_fov_width = math.sqrt(footprint) # FOV_Width
    aircaft_y = ingress_dist_km + (sensor_fov_width/2) # Start in center of fov
    
    while area_sanitized < patrol_area:
        rtb_time = hours_to_fly_xkm(speed_mps, aircaft_y)
        # if debug : print(f"{aircaft_y} km from base ({rtb_time}) hrs" )
        
        #TODO: Account for moving east/west, and changing distance from base
        if endurance <  x_patrol_time+rtb_time: # If we have enough TOS to complete this strip
            print(f"Bingo fuel.  {percent_patrol_complete(patrol_area, area_sanitized):.2f}, % sanitized")
            if debug: print(f"Endurance remaining: {endurance}. Distance from base {aircaft_y}km.  Time for 1 more patrol and RTB:{hours_to_fly_xkm(speed_mps,aoi_width+aircaft_y)} ")
            return area_sanitized
        
        
        endurance                -= x_patrol_time            # deduct "fuel"
        area_sanitized           += sensor_fov_width*aoi_width    # record area sanitized
        aircaft_y     += sensor_fov_width              # move aircraft away from base
        progress_percent          = (area_sanitized / patrol_area) * 100
        # print(f"\rPatrol progress: {progress_percent:.1f}% complete", end="", flush=True)
    assert(hours_to_fly_xkm(speed_mps, aircaft_y) <= endurance)
    print(f"Mission complete. {percent_patrol_complete(patrol_area, area_sanitized):.2f}, % sanitized")
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