from models import *
from simulation.doe import create_doe

import numpy as np

### Factors
# Speed (Mach) from 0.4 to 0.9, e.g., 6 levels
speeds = np.linspace(0.4, 0.9, 6)
# Altitude from 5k to 25k feet, convert to meters (1 ft = 0.3048 m)
altitudes_ft = np.linspace(5000, 25000, 5)
###

doe = create_doe(speeds, altitudes_ft)