import numpy as np

def analytic_ballistic(z0, V0, gamma_deg, g0=9.80665, n_steps=500):
    """
    Returns time, x, z arrays for a vacuum trajectory.
    """
    gamma = np.deg2rad(gamma_deg)
    vx0, vz0 = V0 * np.cos(gamma), V0 * np.sin(gamma)
    # Solve time-of-flight:
    T = ( vz0 + np.sqrt(vz0**2 + 2 * g0 * z0 ) ) / g0
    t = np.linspace(0, T, n_steps)
    x = vx0 * t
    z = z0 + vz0 * t - 0.5 * g0 * t**2
    return t, x, z

# Example usage:
t_ana, x_ana, z_ana = analytic_ballistic(
    z0=70000,      # initial altitude [m]
    V0=2500,       # re-entry speed [m/s]
    gamma_deg=-30  # flight-path angle [deg]
)
