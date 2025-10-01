from __future__ import annotations

import logging
import math
import os
from typing import Tuple

import numpy as np

# Local imports (prefer explicit over wildcard; keep as-is if your module names differ)
# from models import cost_model, endurance_model
from models import *  # noqa: F401,F403  # If you truly need wildcard, leave it and remove the line above.
from simulation.doe import create_doe

# -----------------------------------------------------------------------------
# Config / constants
# -----------------------------------------------------------------------------
METERS_PER_FOOT = 0.3048
KM_PER_METER = 1e-3
G = 9.81  # m/s^2
SEA_LEVEL_TEMP_K = 288.15  # K
LAPSE_RATE_K_PER_M = 0.0065  # K/m
SEA_LEVEL_A_MS = 340.3  # m/s  (speed of sound at 288.15 K)

logging.basicConfig(
    level=logging.INFO,
    format="%(levelname)s:%(name)s:%(message)s",
)
log = logging.getLogger("mpa")


# -----------------------------------------------------------------------------
# Geometry / sensor helpers
# -----------------------------------------------------------------------------
def get_effective_fov_deg(
    alt_ft: float,
    fov_deg: float,
    npx_axis: int,
    target_m: float = 15.0,
    px_req: int = 3,
) -> float:
    """
    Compute effective FOV (deg) so that the pixel requirement holds at the worst-case edge.

    Parameters
    ----------
    alt_ft : float
        Altitude (ft).
    fov_deg : float
        Optical field of view along the axis (deg).
    npx_axis : int
        Pixel count along width.
    target_m : float
        Target dimension along that axis (m). Default 15 (approx. corvette beam).
    px_req : int
        Required pixels across target at the edge. Default 3.

    Returns
    -------
    float
        Effective FOV in degrees (<= fov_deg).
    """
    if alt_ft < 0:
        raise ValueError("Altitude must be non-negative.")
    if npx_axis <= 0:
        raise ValueError("npx_axis must be positive.")
    if target_m <= 0:
        raise ValueError("target_m must be positive.")
    if px_req <= 0:
        raise ValueError("px_req must be positive.")

    h_m = alt_ft * METERS_PER_FOOT
    fov_rad = math.radians(fov_deg)
    ifov_rad_per_px = fov_rad / npx_axis

    # ratio = (px_req * h * IFOV) / target_size  (dimensionless)
    ratio = (px_req * h_m * ifov_rad_per_px) / target_m
    ratio = max(0.0, ratio)  # guard tiny negative from float noise

    if ratio >= 1.0:
        theta_max = 0.0
    else:
        theta_max = math.acos(math.sqrt(ratio))

    fov_eff_rad = min(fov_rad, 2.0 * theta_max)
    return math.degrees(fov_eff_rad)


def get_sensor_footprint_w_km(
    alt_ft: float,
    sensor,
    target_m: float = 15.0,
    px_req: int = 3,
) -> float:
    """
    Compute sensor footprint WIDTH (km) using the *effective* FOV to enforce pixel criteria.

    Parameters
    ----------
    alt_ft : float
        Altitude (ft).
    sensor : object
        Must provide attributes: fov_deg_w (deg), npx_w (px).
    target_m : float
        Target dimension along width axis (m).
    px_req : int
        Required pixels across target at the edge.

    Returns
    -------
    float
        Footprint width in km.
    """
    h_m = alt_ft * METERS_PER_FOOT

    fov_eff_w_deg = get_effective_fov_deg(
        alt_ft=alt_ft,
        fov_deg=sensor.fov_deg_w,
        npx_axis=sensor.npx_w,
        target_m=target_m,
        px_req=px_req,
    )
    width_m = 2.0 * h_m * math.tan(math.radians(fov_eff_w_deg / 2.0 * 2.0))
    # NOTE: math.tan( (FOV_eff/2) * pi/180 ). The *2.0 cancelled because we already split by 2 above.
    # To avoid confusion, write it explicitly:
    width_m = 2.0 * h_m * math.tan(math.radians(fov_eff_w_deg) / 2.0)

    return width_m * KM_PER_METER


# -----------------------------------------------------------------------------
# Flight / performance helpers
# -----------------------------------------------------------------------------
def mach_to_mps(mach: float, alt_ft: float) -> float:
    """
    Compute true airspeed (m/s) from Mach and altitude (ft), using a simple ISA lapse model:

      T(h) = 288.15 - 0.0065 * h(m)
      a(h) = 340.3 * sqrt(T/288.15)

    Parameters
    ----------
    mach : float
        Mach number.
    alt_ft : float
        Altitude (ft).

    Returns
    -------
    float
        Speed (m/s).
    """
    alt_m = alt_ft * METERS_PER_FOOT
    t_kelvin = SEA_LEVEL_TEMP_K - LAPSE_RATE_K_PER_M * alt_m
    # Guard against negative due to model misuse at very high altitudes:
    t_kelvin = max(1.0, t_kelvin)
    a_ms = SEA_LEVEL_A_MS * math.sqrt(t_kelvin / SEA_LEVEL_TEMP_K)
    return mach * a_ms


def hours_to_fly_xkm(speed_mps: float, distance_km: float) -> float:
    """
    Convert distance (km) at speed (m/s) to hours.
    """
    distance_m = distance_km * 1_000.0
    time_sec = distance_m / speed_mps
    return time_sec / 3600.0


def percent_patrol_complete(gridsize_km2: float, area_sanitized_km2: float) -> float:
    return 100.0 * (area_sanitized_km2 / gridsize_km2) if gridsize_km2 > 0 else 0.0


def distance_from_base(aircraft_x_km: float, aircraft_y_km: float, ingress_dist_km: float) -> float:
    """
    Distance (km) from base at (0, -ingress_dist_km) to aircraft position (x, y) in km.
    
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
        """
    dx = aircraft_x_km - 0.0
    dy = aircraft_y_km + ingress_dist_km  # base is SOUTH of origin
    return math.hypot(dx, dy)


def turn_fuel_multiplier(bank_angle_deg: float, k_induced_drag: float = 0.2) -> float:
    """
    Fuel-burn multiplier for a *level* banked turn relative to straight & level.

    multiplier = 1 + k * (n^2 - 1),  with  n = 1 / cos(phi)

    Parameters
    ----------
    bank_angle_deg : float
        Bank angle (deg).
    k_induced_drag : float
        Fraction of total drag at cruise due to induced drag (0..1). Default 0.2.

    Returns
    -------
    float
        Fuel-burn multiplier (>=1).
    """
    phi_rad = math.radians(bank_angle_deg)
    n = 1.0 / max(1e-9, math.cos(phi_rad)) #load factor
    return 1.0 + k_induced_drag * (n ** 2 - 1.0)


def calc_turn_radius(speed_mps: float, bank_angle_deg: float) -> float:
    """
    Turn radius (m) for a coordinated level turn.
    """
    bank_angle_rad = math.radians(bank_angle_deg)
    tan_phi = math.tan(bank_angle_rad)
    if abs(tan_phi) < 1e-12:
        return float("inf")
    return (speed_mps ** 2) / (G * tan_phi)


def calc_orbit_turn_time(
    speed_mps: float,
    sensor_footprint_width_km: float,
    bank_angle_deg: float,
    alt_ft: float,
    speed_mach: float,
) -> Tuple[float, float]:
    """
    Time to accomplish a “lane-change” turn between parallel sanitization tracks.

    If the diameter (2R) is larger than the lane spacing (W), the aircraft must overshoot,
    trace additional arcs, and then rejoin. Otherwise, a simple half-circle suffices.

    Parameters
    ----------
    speed_mps : float
        True airspeed (m/s).
    sensor_footprint_width_km : float
        Lane spacing (km) i.e., effective footprint width W.
    bank_angle_deg : float
        Bank angle (deg).
    alt_ft : float
        For reporting only.
    speed_mach : float
        For reporting only.

    Returns
    -------
    (time_hrs_with_multiplier, time_hrs) : Tuple[float, float]
        Hours with fuel-burn multiplier applied (for endurance accounting),
        and raw geometric time (for timeline).
    """
    radius_m = calc_turn_radius(speed_mps, bank_angle_deg)
    w_m = sensor_footprint_width_km * 1_000.0

    log.debug(
        "Mach: %.3f, alt_ft: %.0f, speed: %.2f m/s, bank: %.1f°, R: %.1f m",
        speed_mach, alt_ft, speed_mps, bank_angle_deg, radius_m,
    )

    if 2.0 * radius_m > w_m:
        # Isosceles triangle with equal sides 2R and base (2R + W)
        base_length = 2.0 * radius_m + w_m
        equal_side = 2.0 * radius_m

        # Vertex angle c via Law of Cosines
        cos_c = (2.0 * equal_side ** 2 - base_length ** 2) / (2.0 * equal_side ** 2)
        cos_c = max(min(cos_c, 1.0), -1.0)
        c = math.acos(cos_c)

        # Base angles
        a = b = (math.pi - c) / 2.0

        # Arcs
        turn_outbound = radius_m * a
        turn_inbound = radius_m * b
        reflex_arc = radius_m * (2.0 * math.pi - c)

        maneuver_distance_m = turn_outbound + turn_inbound + reflex_arc
    else:
        # Simple half-circle turn
        maneuver_distance_m = math.pi * radius_m

    time_hrs = hours_to_fly_xkm(speed_mps, maneuver_distance_m * KM_PER_METER)
    time_hrs_with_mult = time_hrs * turn_fuel_multiplier(bank_angle_deg)
    return time_hrs_with_mult, time_hrs


# -----------------------------------------------------------------------------
# Patrol simulator
# -----------------------------------------------------------------------------
def calc_area_sanitized(
    alt_ft: float,
    speed_mach: float,
    speed_mps: float,
    endurance_hr: float,
    footprint_w_km: float,
    ingress_dist_km: float,
    bank_angle_deg: float = 60.0,
    fuel_reserve_hrs: float = 0.0,
    aoi_width_km: float = 100.0,
    aoi_length_km: float = 100.0,
    verbose: bool = False,
) -> Tuple[float, float]:
    """
    Simulate an MPA sanitizing a rectangular AOI using parallel tracks.
    The aircraft starts at the southern edge of the AOI (0,0) after ingress and flies
    strips northward, deducting endurance for each strip and the return trip to base.

    Base coordinates are at (0, -ingress_dist_km), directly south of AOI origin.

    AOI dimensions:
        Width  = aoi_width  (km, east-west)
        Length = aoi_length (km, north-south)

    Sensor footprint is assumed square

    Coordinates schematic (for ingress_dist_km = 100 km, aoi_width = 100 km):

          y (north)
          ↑
          |
          |       * Aircraft position (aircraft_x, aircraft_y)
          |       * For first strip, aircraft_y = footprint
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
        footprint (float): Sensor footprint width in km.
        ingress_dist_km (float): Distance from base to AOI southern edge (km).
        fuel_reserve_hrs (float): Minimum reserve endurance to retain (hours).
        aoi_width (float, optional): Width of AOI (km, default 100 km).
        aoi_length (float, optional): Length of AOI (km, default 100 km).

    Returns
    -------
    (area_sanitized_km2, mission_time_hrs)
    """
    if bank_angle_deg >= 90.0:
        raise ValueError("Bank angle must be < 90°.")

    # Stopwatch (timeline) and fuel (endurance) bookkeeping
    clock_hrs = 0.0
    endurance_hrs = endurance_hr

    patrol_area_km2 = aoi_width_km * aoi_length_km
    w_km = footprint_w_km

    # Start at west edge, half-footprint inside the AOI
    aircraft_x_km = 0.0
    aircraft_y_km = w_km / 2.0
    x_edges = (0.0, aoi_width_km)

    # One “east-west strip” time (no turn)
    strip_time_hrs = hours_to_fly_xkm(speed_mps, aoi_width_km)

    # Reserve + ingress
    endurance_hrs -= fuel_reserve_hrs
    ingress_time_hrs = hours_to_fly_xkm(speed_mps, ingress_dist_km)
    endurance_hrs -= ingress_time_hrs
    clock_hrs += ingress_time_hrs

    # Precompute per-turn times
    turn_time_hrs_fuel, turn_time_hrs = calc_orbit_turn_time(
        speed_mps, w_km, bank_angle_deg, alt_ft, speed_mach
    )

    area_sanitized_km2 = 0.0

    if verbose:
        log.info("Starting endurance (after reserve + ingress): %.3f h", endurance_hrs)

    while area_sanitized_km2 < patrol_area_km2:
        # Time to RTB from current location
        km_from_base = distance_from_base(aircraft_x_km, aircraft_y_km, ingress_dist_km)
        rtb_time_hrs = hours_to_fly_xkm(speed_mps, km_from_base)

        # If we cannot finish the next strip AND get home, RTB now
        if endurance_hrs < strip_time_hrs + rtb_time_hrs:
            if verbose:
                log.info(
                    "Bingo fuel. %.2f%% sanitized.",
                    percent_patrol_complete(patrol_area_km2, area_sanitized_km2),
                )
            clock_hrs += rtb_time_hrs
            return area_sanitized_km2, clock_hrs

        # Fly the strip + turn for next lane
        endurance_hrs -= strip_time_hrs
        endurance_hrs -= turn_time_hrs_fuel

        area_sanitized_km2 += w_km * aoi_width_km
        clock_hrs += strip_time_hrs + turn_time_hrs

        # Move to next lane (toggle east/west, step north by W)
        aircraft_x_km = x_edges[1] if aircraft_x_km == x_edges[0] else x_edges[0]
        aircraft_y_km += w_km

        if verbose:
            log.debug(
                "Pos=(%.1f,%.1f) km, sanitized=%.1f km², endurance=%.2f h",
                aircraft_x_km, aircraft_y_km, area_sanitized_km2, endurance_hrs,
            )

        # If we exit the AOI (fully covered), break to RTB
        if aircraft_y_km >= aoi_length_km + w_km / 2.0:
            break

    # AOI complete: compute RTB from final position
    km_from_base_final = distance_from_base(aircraft_x_km, aircraft_y_km, ingress_dist_km)
    rtb_time_hrs_final = hours_to_fly_xkm(speed_mps, km_from_base_final)

    # If this fails, the AOI was “complete” but fuel won’t get you home; caller should handle policy.
    # if endurance_hrs < rtb_time_hrs_final:
    #     log.warning(
    #         "AOI complete, but insufficient endurance for RTB: need %.2f h, have %.2f h",
    #         rtb_time_hrs_final, endurance_hrs,
    #     )

    clock_hrs += rtb_time_hrs_final
    return patrol_area_km2, clock_hrs


# -----------------------------------------------------------------------------
# Main DOE driver
# -----------------------------------------------------------------------------
def main() -> None:
    # Assumptions
    ingress_dist_km = 100.0
    fuel_reserve_hrs = 1.0
    bank_angle_deg = 45.0

    # Factors
    machs = np.linspace(0.4, 0.9, 6)
    alts_ft = np.linspace(5_000.0, 25_000.0, 5)

    doe = create_doe(machs, alts_ft)

    # Derived columns
    doe["speed_mps"] = doe.apply(lambda r: mach_to_mps(r.speed_mach, r.alt_ft), axis=1)
    doe["cost"] = doe.apply(lambda r: cost_model(r.speed_mach, r.alt_ft / 1_000.0, r.sensor), axis=1)
    doe["endurance_hr"] = doe.apply(lambda r: endurance_model(r.speed_mach, r.alt_ft / 1_000.0), axis=1)
    doe["footprint_w_km"] = doe.apply(lambda r: get_sensor_footprint_w_km(r.alt_ft, r.sensor), axis=1)

    results_area = []
    results_time = []

    for _, row in doe.iterrows():
        area_km2, time_hrs = calc_area_sanitized(
            alt_ft=row["alt_ft"],
            speed_mach=row["speed_mach"],
            speed_mps=row["speed_mps"],
            endurance_hr=row["endurance_hr"],
            footprint_w_km=row["footprint_w_km"],
            ingress_dist_km=ingress_dist_km,
            bank_angle_deg=bank_angle_deg,
            fuel_reserve_hrs=fuel_reserve_hrs,
            aoi_width_km=100.0,
            aoi_length_km=100.0,
            verbose=False,
        )
        results_area.append(area_km2)
        results_time.append(time_hrs)

    doe["area_sanitized_km2"] = results_area
    doe["mission_time_hrs"] = results_time

    # Ensure output dir exists
    os.makedirs("output", exist_ok=True)
    # Do not write the raw sensor object to CSV
    cols = [c for c in doe.columns if c != "sensor"]
    doe.to_csv("output/mpa.csv", columns=cols, index=False)


if __name__ == "__main__":
    main()
