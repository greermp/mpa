from models import *
from simulation.doe import create_doe

import numpy as np

def main():

    ### Factors
    # Speed (Mach) from 0.4 to 0.9, e.g., 6 levels
    machs = np.linspace(0.4, 0.9, 6)
    # Altitude from 5k to 25k feet
    altitudes_ft = np.linspace(5000, 25000, 5)
    ###
    
    # Assume a spacing between sweeps of footprint * (1-sweep_overlap) 
    sweep_overlap = 0.2

    doe = create_doe(machs, altitudes_ft)
    doe['cost'] = doe.apply(
        lambda row: cost_model(row.speed_mach, row.alt_ft, row.sensor), axis=1
    )
    doe['endurance_hr'] = doe.apply(
        lambda row: endurance_model(row.speed_mach, row.alt_ft), axis=1
    )





if __name__ == "__main__":
    main()