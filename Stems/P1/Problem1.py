import numpy as np
from scipy.integrate import quad

def calculate_sto(thrust, weight, cl_max, cd, rho, s, gc):
    def v_stall(weight, rho, s, cl_max):
        return np.sqrt(weight / (0.5 * rho * s * cl_max))

    def sto_integral(v, a, b):
        return v / (a - b * v**2)

    v_stall_value = v_stall(weight, rho, s, cl_max)
    v_to = 1.2 * v_stall_value
    a = (gc * thrust) / weight
    b = (gc / weight) * (0.5 * rho * s * cd)

    sto_distance, _ = quad(sto_integral, 0, v_to, args=(a, b))
    return sto_distance

# Example parameters
rho = 0.002377  # slugs/ft^3 (density at sea level)
s = 1000  # wing area in ft^2
cl_max = 2.4
cd = 0.0279
gc = 32.2  # ft/s^2

# Example usage
thrust = 13000  # in lbs
weight = 56000  # in lbs
sto_distance = calculate_sto(thrust, weight, cl_max, cd, rho, s, gc)
print(f"Calculated STO: {sto_distance} feet")
