from models import *
from simulation.doe import create_doe

import math
import sys

import numpy as np
debug = False

def get_sensor_footprint_km2(alt_ft, sensor):
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
        Footprint in km2.
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
    T0 = 288.15  # sea level temperature in K
    lapse_rate = 0.0065
    alt_m    = alt_ft * 0.3048
    t_kelvin = T0 - lapse_rate * alt_m
    a        = 340.3 * math.sqrt(t_kelvin/T0) # a: Speed of Sound (m/s)
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

def distance_from_base(aircraft_x, aircraft_y, ingress_dist, side = "W"):
    """
    Compute distance to base (0, -ingress_dist) to aircraft position
                y (north)
                ↑
                |
                
             _ _ _ _ _ _ 
            |           |
            |           |
            |           |  
AOI origin  |           |
(0,0)       *-----------* <- first patrol leg
            W           E
            |
            *
             Base (0, -ingress_dist)

    Parameters
    ----------
    aircraft_x : float
        x-position within AOI (km, east-west)
    aircraft_y : float
        y-position within AOI (km, north)
    ingress_dist : float
        Distance from base to AOI southern edge (km)
    side : str, optional
        Side of the AOI the patrol starts on, "W" or "E" (default "W")

    Returns
    -------
    float
        Distance from base to aircraft (km)
    """
    base_x = 0
    base_y = -ingress_dist  # base south of AOI origin
    dx = aircraft_x - base_x
    dy = aircraft_y - base_y
    return math.hypot(dx, dy)

def turn_fuel_multiplier(bank_angle_deg, k_induced_drag=0.2):
    """
    Calculate the fuel burn multiplier for a level banked turn.

    Parameters:
        bank_angle_deg (float): Bank angle in degrees.
        k_induced_drag (float): Fraction of total drag that is induced drag at cruise (default 0.2).

    Returns:
        float: Fuel burn multiplier relative to straight & level flight.
    """
    # Convert bank angle to radians
    phi_rad = math.radians(bank_angle_deg)

    # Load factor n
    n = 1 / math.cos(phi_rad)

    # Fuel burn multiplier
    multiplier = 1 + k_induced_drag * (n**2 - 1)

    return multiplier


def calc_orbit_turn_time(speed_mps, sensor_footprint_width_km, bank_angle_deg, alt_ft, speed_mach):
    """
    Compute the time to fly a turn maneuver given aircraft speed, sensor footprint width, and bank angle.
    Since we are utilzing endurance as a surrogate for fuel, also return time*scalar to account for more
    drag and thrust required (Fuel flow)

    Scenario:
    - If the aircraft’s turn radius is large compared to the sensor footprint width,
      the aircraft must overshoot the turn, fly an arc around the turn circle center,
      and then rejoin the patrol leg.
    - Otherwise, a simple half-circle turn is sufficient.
    
                     X
                                          
                     P    
                    / \    
                  /     \   
         b      /         \     c
       .      /             \
      .      /               \
    R.      /                 \
    . c     *----x--------y-----* <- first patrol leg
   ..........        |            
        R            |
            Lane X   | Lane Y
                     |  
            *
             Base (0, -ingress_dist)

    Parameters
    ----------
    speed_mps : float
        True airspeed of the aircraft in meters/second.
    sensor_footprint_width_km : float
        Sensor footprint width (meters).
    bank_angle_deg : float
        Bank angle of the aircraft (degrees).


    Returns
    -------
    time_to_turn_hrs_wmult : float
        Time in seconds to fly from x to y
    time_to_turn : float
        Time in seconds to fly from x to y
    """
    radius_meters = calc_turn_radius(speed_mps, bank_angle_deg)
    sensor_footprint_width_meters = sensor_footprint_width_km * 1000
    
    print(f"Mach: {speed_mach}, alt: {alt_ft}, mps: {speed_mps}, bank: {bank_angle_deg}, turn radius {radius_meters}")
    # # Can we slow down??
    # if (2 * radius_meters > sensor_footprint_width_meters):
    #     speed_mps = mach_to_mps(0.5, alt_ft)
    #     radius_meters = calc_turn_radius(speed_mps, bank_angle_deg)
    #     sensor_footprint_width_meters = sensor_footprint_width_km * 1000
    
    if (2 * radius_meters > sensor_footprint_width_meters):
        # print(f"Turn diameter {2*radius_meters} > sensor footprint width {sensor_footprint_width_meters}")

        base_length = 2 * radius_meters + sensor_footprint_width_meters
        equal_side = 2 * radius_meters

        # Vertex angle c using Law of Cosines
        cos_c = (2 * equal_side**2 - base_length**2) / (2 * equal_side**2)
        cos_c = max(min(cos_c, 1), -1)
        c = math.acos(cos_c)

        # Base angles a and b (equal in isosceles triangle)
        a = b = (math.pi - c) / 2

        # print(f"Base angles: a = {math.degrees(a):.2f}°, b = {math.degrees(b):.2f}°, vertex angle c = {math.degrees(c):.2f}°")

        # Arc lengths for outbound and inbound segments
        turn_outbound = radius_meters * a
        turn_inbound = radius_meters * b  # a == b

        # print(f"Outbound arc length: {turn_outbound:.2f} m")
        # print(f"Inbound arc length: {turn_inbound:.2f} m")

        # Reflex angle and corresponding arc for the main turn
        reflex_angle = 2 * math.pi - c
        reflex_arc = radius_meters * reflex_angle

        # print(f"Reflex angle: {math.degrees(reflex_angle):.2f}°, arc length: {reflex_arc:.2f} m")

        # Total maneuver distance
        maneuver_distance = turn_outbound + turn_inbound + reflex_arc
        print(f"Total maneuver distance: {maneuver_distance:.2f} m")
    else:
        # Simple half-circle turn
        maneuver_distance = math.pi * radius_meters
        print(f"Simple half-circle turn: distance = {maneuver_distance:.2f} m")

    # Time = distance / speed
    time_to_turn_hrs = hours_to_fly_xkm(speed_mps, maneuver_distance/1000)
    time_to_turn_hrs_wmult = time_to_turn_hrs * turn_fuel_multiplier(bank_angle_deg)
    
    # print(f"Time to turn: {maneuver_distance/1000:.2f} km  {time_to_turn_hrs:.2f} h. With multiplier {calculate_turn_fuel_penalty(bank_angle_deg, 'n15')}: {time_to_turn_hrs_wmult}hrs")
    # print(f"Time to patrol 100km: {hours_to_fly_xkm(speed_mps, 100)} h")
    return time_to_turn_hrs_wmult, time_to_turn_hrs

def calc_turn_radius(speed_mps, bank_angle_deg):
    """
    Compute turn radius (m) given speed and bank angle.
    """
    g = 9.81  # m/s^2
    bank_angle_rad = math.radians(bank_angle_deg)
    return speed_mps**2 / (g * math.tan(bank_angle_rad))
    
def calc_area_sanitized(
    alt_ft,
    speed_mach,
    speed_mps,
    endurance,
    sensor_name,
    footprint,
    ingress_dist_km,
    # Assumptions
    bank_angle_deg=60,
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
          |       * For first strip, aircraft_y = math.sqrt(footprint)/2
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
        float: Time to complete patrol
    """
    if bank_angle_deg >= 90:
        raise ValueError("Bank angle must be less than 90°")
    
    # Mission
    stop_watch = 0
    patrol_area    = aoi_width*aoi_length
    sensor_footprint_width = math.sqrt(footprint) # FOV_Width

    # Geometry    
    aircraft_y     = sensor_footprint_width/2 # Start laterally displaced half sensor fov (optimal)
    aircraft_x     = 0 # 0 representes W side of the grid, 1 the E side
    x_edges        = [0, aoi_width]  # west and east edges
    
    # MOP
    area_sanitized = 0
    x_patrol_time  = hours_to_fly_xkm(speed_mps, aoi_width) # time to sanitize from position (0, n) to (aoi_width, n)
    
    # Fuel accouting
    endurance     -= fuel_reserve_hrs # Subtract reserve fuel from `endurance`
    endurance     -= hours_to_fly_xkm(speed_mps, ingress_dist_km) # Subtract time to ingress from `endurance`
    stop_watch    += hours_to_fly_xkm(speed_mps, ingress_dist_km)
    
    # Calculate time to turn. Function of turn radius and sensor footprint
    # For fuel  ........... for time
    time_to_turn_hrs_wmult, time_to_turn_hrs = calc_orbit_turn_time(speed_mps, sensor_footprint_width, bank_angle_deg, alt_ft, speed_mach)

    print (f"Starting Endurance : {endurance}")
    while area_sanitized < patrol_area:
        km_from_base = distance_from_base(aircraft_x, aircraft_y, ingress_dist_km)
        if debug: print(f"X:{aircraft_x}, Y:{aircraft_y},  Distance to base:{km_from_base}")
        rtb_time = hours_to_fly_xkm(speed_mps, km_from_base)
        # if debug : print(f"{aircraft_y} km from base ({rtb_time}) hrs" )
        
        if endurance <  x_patrol_time+rtb_time: # If we don't have enough TOS to complete this strip
            print(f"Bingo fuel.  {percent_patrol_complete(patrol_area, area_sanitized):.2f}, % sanitized")
            if debug: print(f"Endurance remaining: {endurance}. Distance from base {aircraft_y}km.  Time for 1 more patrol and RTB:{hours_to_fly_xkm(speed_mps,aoi_width+aircraft_y)} ")
            return area_sanitized, (stop_watch+rtb_time)

        
        
        print (f"Deduct patrol time : {x_patrol_time}")
        endurance                -= x_patrol_time 
        print (f"Deduct turn time : {time_to_turn_hrs_wmult}")
        endurance                -= time_to_turn_hrs_wmult                 # deduct "fuel"
        print (f"Endurance : {endurance}")
        area_sanitized           += sensor_footprint_width * aoi_width    # record area sanitized
        stop_watch               += time_to_turn_hrs + x_patrol_time
        aircraft_y               += sensor_footprint_width              # move aircraft away from base
        aircraft_x = x_edges[1] if aircraft_x == x_edges[0] else x_edges[0] # Track W (0) or East (100)
        print(f"Aircraft at ({aircraft_x},{aircraft_y}).  Sanitized: {area_sanitized}")
    assert(hours_to_fly_xkm(speed_mps, aircraft_y) <= endurance)
    print(f"Mission complete. {percent_patrol_complete(patrol_area, area_sanitized):.2f}, % sanitized")
    return patrol_area, (stop_watch+rtb_time)
            
def main():
    ### Assumptions
    ingress_dist_km  = 100
    fuel_reserve_hrs =   1
    bank_angle_deg   =  45
    
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
        lambda row: get_sensor_footprint_km2(row.alt_ft, row.sensor), axis=1)

    for idx, row in doe.iterrows():
        area_sanitized_km2, time_hrs = calc_area_sanitized(row["alt_ft"],
                                                                 row["speed_mach"],
                                                                 row["speed_mps"],
                                                                 row["endurance_hr"],
                                                                 row["sensor_name"],
                                                                 row["footprint_km2"],
                                                                 ingress_dist_km=ingress_dist_km, 
                                                                 bank_angle_deg = bank_angle_deg,
                                                                 fuel_reserve_hrs=fuel_reserve_hrs,
                                                                 aoi_width=100,
                                                                 aoi_length=375
                                                                 )
        doe.loc[idx, "area_sanitized_km2"] = area_sanitized_km2
        doe.loc[idx, "mission_time_hrs"] = time_hrs
        # break
    
    doe.to_csv("output/mpa_big.csv", columns= [col for col in doe.columns if col != "sensor"] ,index=False)

if __name__ == "__main__":
    main()
    # mach   = 0.6
    # mps   = mach_to_mps(mach, 10000)
    # calc_orbit_turn_time(mps, 4000, 45)
    # print(f"Mach: {mach} mps: {mps} Radius: {radius} Turn time: {turn_time(mps, radius)*60} minutes")