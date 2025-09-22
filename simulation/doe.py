from models.sensor import sensors

from itertools import product

import pandas as pd

def create_doe(machs, altitudes_ft):
    doe = pd.DataFrame(list(product(machs, altitudes_ft, sensors)), 
                       columns = ['speed_mach', 'alt_ft', 'sensor'])
    doe['sensor_name'] = doe['sensor'].apply(lambda s: s.name)
    return doe