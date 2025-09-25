import pytest
from models import *
import math
import pytest
from sim import get_sensor_footprint

@pytest.fixture
def sensor15():
    """Returns a lightweight 15° x 15° sensor"""
    return EO_Sensor.from_fov(15, 15,0.05)

@pytest.fixture
def sensor30():
    """Returns a lightweight 30° x 30° sensor"""
    return EO_Sensor.from_fov(30, 30, 1)

@pytest.fixture
def sensor60():
    """Returns a lightweight 60° x 60° sensor"""
    return EO_Sensor.from_fov(60, 60, 10)

def test_15deg_5000ft(sensor15):
    result = get_sensor_footprint(5000, sensor15)
    assert math.isclose(result, 0.161, rel_tol=1e-2)

def test_60deg_25000ft(sensor60):
    result = get_sensor_footprint(25000, sensor60)
    assert math.isclose(result, 77.42, rel_tol=1e-2)
    
    
# -----------------------------
# Sensor footprint tests
# -----------------------------
@pytest.mark.parametrize("alt_ft,sensor_fixture,expected", [
    (5000, "sensor15", 0.161),
    (10000, "sensor30", 2.66),
    (25000, "sensor60", 77.42),
])
def test_sensor_footprint(alt_ft, sensor_fixture, expected, request):
    sensor = request.getfixturevalue(sensor_fixture)
    result = get_sensor_footprint(alt_ft, sensor)
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