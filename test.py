import pytest
import math

from models import *
from sim import (
    get_sensor_footprint_km2,
    cost_model,
    endurance_model,
    calc_orbit_turn_time,
    calc_turn_radius,
    mach_to_mps,
    calculate_turn_fuel_penalty,
    hours_to_fly_xkm
)

# -----------------------------
# Fixtures
# -----------------------------
@pytest.fixture
def sensor15():
    """Returns a lightweight 15° x 15° sensor"""
    return EO_Sensor.from_fov(15, 15, 0.05)

@pytest.fixture
def sensor30():
    """Returns a lightweight 30° x 30° sensor"""
    return EO_Sensor.from_fov(30, 30, 1)

@pytest.fixture
def sensor60():
    """Returns a lightweight 60° x 60° sensor"""
    return EO_Sensor.from_fov(60, 60, 10)


# -----------------------------
# Sensor footprint tests
# -----------------------------
def test_15deg_5000ft(sensor15):
    result = get_sensor_footprint_km2(5000, sensor15)
    assert math.isclose(result, 0.161, rel_tol=1e-2)

def test_60deg_25000ft(sensor60):
    result = get_sensor_footprint_km2(25000, sensor60)
    assert math.isclose(result, 77.42, rel_tol=1e-2)

@pytest.mark.parametrize("alt_ft,sensor_fixture,expected", [
    (5000, "sensor15", 0.161),
    (10000, "sensor30", 2.66),
    (25000, "sensor60", 77.42),
])
def test_sensor_footprint(alt_ft, sensor_fixture, expected, request):
    sensor = request.getfixturevalue(sensor_fixture)
    result = get_sensor_footprint_km2(alt_ft, sensor)
    assert math.isclose(result, expected, rel_tol=1e-2)


# -----------------------------
# Cost model tests
# -----------------------------
@pytest.mark.parametrize("mach, alt_kft, sensor_fixture, expected", [
    (0.5, 10, "sensor30", 50*0.5**2 - 35*0.5 + 0.03*10**2 - 0.02*10 + 11 + 1),
    (0.4, 5,  "sensor15", 50*0.4**2 - 35*0.4 + 0.03*5**2  - 0.02*5  + 11 + 0.05),
    (0.8, 25, "sensor60", 50*0.8**2 - 35*0.8 + 0.03*25**2 - 0.02*25 + 11 + 10.0)
])
def test_cost_model(mach, alt_kft, sensor_fixture, expected, request):
    sensor = request.getfixturevalue(sensor_fixture)
    result = cost_model(mach, alt_kft, sensor)
    assert math.isclose(result, expected, rel_tol=1e-5)


# -----------------------------
# Endurance model tests
# -----------------------------
@pytest.mark.parametrize("mach, alt_kft, expected", [
    (0.4, 5,  -18.75*0.4**2 + 8.0893*0.4 + 0.01*5**2  + 0.05*5  + 9.2105),
    (0.5, 10, -18.75*0.5**2 + 8.0893*0.5 + 0.01*10**2 + 0.05*10 + 9.2105),
    (0.8, 25, -18.75*0.8**2 + 8.0893*0.8 + 0.01*25**2 + 0.05*25 + 9.2105)
])
def test_endurance_model(mach, alt_kft, expected):
    result = endurance_model(mach, alt_kft)
    assert math.isclose(result, expected, rel_tol=1e-5)


# -----------------------------
# Turn time model tests
# -----------------------------
@pytest.mark.parametrize("speed_mps,sensor_width_km,bank_angle_deg,alt_ft,speed_mach", [
    (200, 1.8, 45, 5000, 0.6),   # Large turn radius scenario
    (150, 0.5, 30, 10000, 0.4),  # Small turn radius scenario
    (250, 2.5, 60, 25000, 0.8),  # Very large turn radius
])
def test_calc_orbit_turn_time(speed_mps, sensor_width_km, bank_angle_deg, alt_ft, speed_mach):
    time_result = calc_orbit_turn_time(speed_mps, sensor_width_km, bank_angle_deg, alt_ft, speed_mach)
    assert isinstance(time_result, float)
    assert time_result > 0

    time_steep = calc_orbit_turn_time(speed_mps, sensor_width_km, 60, alt_ft, speed_mach)
    time_shallow = calc_orbit_turn_time(speed_mps, sensor_width_km, 30, alt_ft, speed_mach)
    assert time_steep <= time_shallow + 1e-3


@pytest.mark.parametrize("speed_mps,bank_angle_deg", [
    (200, 30),
    (150, 45),
    (250, 60),
])
def test_calc_turn_radius_consistency(speed_mps, bank_angle_deg):
    radius = calc_turn_radius(speed_mps, bank_angle_deg)
    assert radius > 0
    assert isinstance(radius, float)


def test_turn_time_sensor_width_dependency():
    speed = 200
    bank = 45
    alt = 5000
    mach = 0.6

    time_narrow = calc_orbit_turn_time(speed, 0.5, bank, alt, mach)
    time_wide = calc_orbit_turn_time(speed, 2.0, bank, alt, mach)

    assert time_wide <= time_narrow + 1e-3
