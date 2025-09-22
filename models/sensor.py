import math
class EO_Sensor:
    def __init__(self, name: str, cost: float, fov_deg_w: float, fov_deg_l: float, range_max: float):
        if not isinstance(cost, (int, float)):
            raise TypeError(f"cost must be a number, got {type(cost).__name__}")
        if cost < 0:
            raise ValueError("fixed_cost cannot be negative")
        
        self.name      = name
        self.cost      = cost
        self.fov_deg_w = fov_deg_w
        self.fov_deg_l = fov_deg_l
        self.range_max = range_max
        
    @classmethod
    def from_fov(cls, fov_deg_w, fov_deg_l, cost):
        """
        Lightweight constructor for testing.
        Defaults cost/weight/power to zero.
        """
        return cls(name="", cost=cost, fov_deg_w=fov_deg_w, fov_deg_l=fov_deg_l, range_max=0)
        
sensors = [
    EO_Sensor("EO/IR Sensor 1", 50000/1e6, 15, 15, math.inf),
    EO_Sensor("EO/IR Sensor 2", 1e6/1e6,   30, 30, math.inf),
    EO_Sensor("EO/IR Sensor 3", 10e6/1e6,  60, 60, math.inf)
]