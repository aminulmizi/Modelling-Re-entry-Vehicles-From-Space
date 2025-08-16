import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import argparse
from scipy.interpolate import interp1d

# =============================================================================
# IN VACUO TEST - TOGGLE DO_INVACUO_TEST TO 'TRUE', ENSURE DRAG IS SET TO 0
# =============================================================================
DO_INVACUO_TEST = False

def analytic_ballistic(z0, V0, gamma_deg, g0=9.80665, n_steps=None, t_numeric=None):
    """
    Computes the pure‑vacuum trajectory under constant gravity.
    If t_numeric is given, returns values at exactly those times; otherwise,
    returns n_steps points between t=0 and impact.
    """
    gamma = math.radians(gamma_deg)
    vx0, vz0 = V0 * math.cos(gamma), V0 * math.sin(gamma)
    # time‑of‑flight until z=0:
    T = ( vz0 + math.sqrt(vz0**2 + 2 * g0 * z0) ) / g0

    if t_numeric is not None:
        t = t_numeric
    else:
        t = np.linspace(0, T, n_steps if n_steps else 500)

    x = vx0 * t
    z = z0 + vz0 * t - 0.5 * g0 * t**2
    return t, x, z


# =============================================================================
# CONFIGURATION / PARAMETERS
# =============================================================================
# These are the default parameters. You can override them using command-line args.
DT = 1                          # Time step in seconds
T_FINAL = 300                   # Maximum simulation time in seconds

# Vehicle parameters
MASS = 25000                    # Mass in kg
DIAMETER = 3.66                 # Diameter in meters
AREA = math.pi * (DIAMETER / 2)**2  # Reference area in m²

# Initial conditions (altitude in meters, speed in m/s, flight path angle in degrees)
INITIAL_ALTITUDE = 70000        # in meters
REENTRY_SPEED = 2500            # in m/s
FLIGHT_PATH_ANGLE_DEG = -30     # negative for descent

# Aero data CSV file (should be in the same folder)
AERO_CSV_FILE = 'Falcon9_aero_data1.csv'
# Optionally, output the trajectory data to a file
OUTPUT_TRAJECTORY_CSV = 'trajectory_data.csv'

SHOW_DEBUG_PRINTS = False
# =============================================================================
# DEBUG PRINT FUNCTION
# =============================================================================
def debug_print(*args):
    if SHOW_DEBUG_PRINTS:
        print("DEBUG:", *args)

# =============================================================================
# ATMOSPHERIC PROPERTIES
# =============================================================================
def atmospheric_properties(altitude):
    """
    Calculate atmospheric properties at a given altitude.
    Returns:
        rho: Density [kg/m^3]
        a: Speed of sound [m/s]
    """
    R_specific = 0.287053  # [kJ/(kg·K)]
    gamma = 1.4

    if altitude <= 11000:
        T = 288.15 - 0.0065 * altitude
        P = 101.325 * ((T / 288.15) ** 5.256)
        debug_print(f"Troposphere: alt={altitude}, T={T:.2f}K, P={P:.2f}kPa")
    elif altitude <= 25000:
        T = 216.65
        P = 22.65 * math.exp(1.73 - 0.000157 * altitude)
        debug_print(f"Lower Stratosphere: alt={altitude}, T={T:.2f}K, P={P:.2f}kPa")
    else:
        T = 141.94 + 0.00299 * altitude
        P = 2.488 * ((T / 216.65) ** -11.388)
        debug_print(f"Upper Atmosphere: alt={altitude}, T={T:.2f}K, P={P:.2f}kPa")
    
    rho = P / (R_specific * T)
    a = math.sqrt(gamma * 1000 * R_specific * T)
    debug_print(f"At alt={altitude}: rho={rho:.4f}, a={a:.2f} m/s")
    return rho, a
# =============================================================================
# GRAVITY FUNCTION
# =============================================================================
def gravity_at_altitude(altitude):
    g0 = 9.80665 #m/s^2
    Er = 6371000 # Earths Mean Radious in meters
    return g0 * (Er / (Er + altitude))**2
# =============================================================================
# PIECEWISE AERO PROPERTIES LOOKUP
# =============================================================================
def get_piecewise_aero_properties(altitude, aero_df):
    """
    Return the piecewise constant aero properties (MACH, AoA, CD) based on the current altitude.
    The CSV rows are milestones (e.g., 70000, 55000, etc.), and the values remain in effect
    until the vehicle reaches the next milestone.
    """
    aero_sorted = aero_df.sort_values(by='ALT', ascending=False).reset_index(drop=True)
    debug_print("Aero milestones (sorted):", aero_sorted['ALT'].values)
    
    chosen_row = aero_sorted.iloc[0]
    for idx, row in aero_sorted.iterrows():
        if altitude <= row['ALT']:
            chosen_row = row
            debug_print(f"At alt={altitude}, using milestone {row['ALT']}: MACH={row['MACH']}, AoA={row['AoA']}, CD={row['CD']}")
        else:
            break
    return float(chosen_row['MACH']), float(chosen_row['AoA']), float(chosen_row['CD'])

def get_drag_coefficient(velocity, altitude, aero_df):
    """
    Return the drag coefficient (Cd) using a piecewise lookup from the aero data.
    """
    mach, aoa, Cd = get_piecewise_aero_properties(altitude, aero_df)
    debug_print(f"At alt={altitude}, piecewise values: MACH={mach}, AoA={aoa}, Cd={Cd}")
    return Cd
# =============================================================================
# INTERPOLATION FUNCTION
# =============================================================================
def get_interpolated_drag_coefficient(altitude, aero_df):
    """
    Return the drag coefficient (Cd) using linear interpolation from the aero data.
    """
    # Sort data by altitude
    aero_sorted = aero_df.sort_values(by='ALT')
    altitudes = aero_sorted['ALT'].values
    cds = aero_sorted['CD'].values

    # Create the interpolator
    cd_interp = interp1d(altitudes, cds, kind='linear', fill_value='extrapolate')
    Cd = float(cd_interp(altitude))
    debug_print(f"At alt={altitude:.1f} m, interpolated Cd={Cd:.4f}")
    return Cd
# =============================================================================
# RUNGE-KUTTA 4th ORDER STEP
# =============================================================================
def rk4_step(state, dt, derivs, *args):
    k1 = derivs(state, *args)
    debug_print("k1:", k1)
    k2 = derivs(state + 0.5 * dt * k1, *args)
    debug_print("k2:", k2)
    k3 = derivs(state + 0.5 * dt * k2, *args)
    debug_print("k3:", k3)
    k4 = derivs(state + dt * k3, *args)
    debug_print("k4:", k4)
    state_next = state + dt * (k1 + 2*k2 + 2*k3 + k4) / 6
    debug_print("New state after RK4 step:", state_next)
    return state_next

# =============================================================================
# DERIVATIVES (3DOF)
# =============================================================================
def derivs(state, mass, area, aero_df):
    x, y, z, vx, vy, vz = state
    velocity = math.sqrt(vx**2 + vy**2 + vz**2)
    debug_print(f"At state={state}, speed={velocity:.2f} m/s")

    g = gravity_at_altitude(z)
    gravity = np.array([0, 0, -g])
    
    rho, _ = atmospheric_properties(z)
    debug_print(f"At alt={z}: rho={rho:.4f}")
    
    Cd = get_interpolated_drag_coefficient(z, aero_df)
    debug_print(f"Cd={Cd}")
    
    F_drag = 0.5 * rho * velocity**2 * Cd * area
    debug_print(f"F_drag={F_drag:.2f} N")

    a_drag = np.array([0, 0, 0])
    if velocity > 0:
        a_drag = -F_drag / mass * np.array([vx, vy, vz]) / velocity
    acceleration = gravity + a_drag
    debug_print(f"a_drag={a_drag}, acceleration={acceleration}")

    return np.array([vx, vy, vz, acceleration[0], acceleration[1], acceleration[2]])

# =============================================================================
# SIMULATION
# =============================================================================
def simulate_trajectory(initial_state, mass, area, aero_df, dt, t_final):
    t = 0.0
    state = initial_state.copy()
    t_arr = [t]
    state_arr = [state.copy()]

    while state[2] > 0 and t < t_final:
        speed = math.sqrt(state[3]**2 + state[4]**2 + state[5]**2)
        debug_print(f"Time={t:.1f}s, state={state}, speed={speed:.2f} m/s")
        state = rk4_step(state, dt, derivs, mass, area, aero_df)
        t += dt

        if np.isnan(state).any():
            print(f"ERROR: NaN found at time={t:.1f}s.")
            break

        t_arr.append(t)
        state_arr.append(state.copy())

    return np.array(t_arr), np.array(state_arr)

def interpolate_to_zero_altitude(t_arr, state_arr):
    # Find the last index where altitude is positive
    for i in range(len(state_arr)-2, -1, -1):
        if state_arr[i][2] > 0 and state_arr[i+1][2] <= 0:
            break
    t1, t2 = t_arr[i], t_arr[i+1]
    s1, s2 = state_arr[i], state_arr[i+1]
    z1, z2 = s1[2], s2[2]
    alpha = z1 / (z1 - z2)
    t0 = t1 + alpha * (t2 - t1)
    state0 = s1 + alpha * (s2 - s1)
    return t0, state0, i
# =============================================================================
# MAIN SETUP
# =============================================================================
def main():
    global DEBUG  # Allow modification from CLI arguments if needed
    
    # Load aero data from CSV file
    aero_data = pd.read_csv(AERO_CSV_FILE)
    debug_print("Loaded CSV data:\n", aero_data)
    
    # Convert flight path angle to radians
    flight_path_angle = math.radians(FLIGHT_PATH_ANGLE_DEG)
    
    # Calculate initial velocities
    vx0 = REENTRY_SPEED * math.cos(flight_path_angle)
    vz0 = REENTRY_SPEED * math.sin(flight_path_angle)
    
    # Set initial state: [x, y, altitude, vx, vy, vz]
    initial_state = np.array([0.0, 0.0, INITIAL_ALTITUDE, vx0, 0.0, vz0])
    debug_print("Initial state:", initial_state)
    debug_print(f"Mass={MASS}, Diameter={DIAMETER}, Area={AREA:.2f}")
    
    # Run the simulation
    time_arr, state_arr = simulate_trajectory(initial_state, MASS, AREA, aero_data, DT, T_FINAL)
    
# INVACUO ANALYSIS - ANALYTICAL
    if DO_INVACUO_TEST:
    # Extract numerical trajectory
        x_num = state_arr[:, 0]
        z_num = state_arr[:, 2]

        # === Analytic vacuum solution ===
        t_ana, x_ana, z_ana = analytic_ballistic(
        z0=INITIAL_ALTITUDE,
        V0=REENTRY_SPEED,
        gamma_deg=FLIGHT_PATH_ANGLE_DEG,
        t_numeric=time_arr
        )

        # === Error calculation ===
        dx = x_num - x_ana
        dz = z_num - z_ana

        print("\n====== Trajectory Comparison (Vacuum) ======")
        print(f"Max downrange error: {np.max(np.abs(dx)):.3f} m")
        print(f"Max vertical   error: {np.max(np.abs(dz)):.3f} m")

        # === Time of flight comparison ===
        T_num = time_arr[-1]
        g0 = 9.80665
        gamma = math.radians(FLIGHT_PATH_ANGLE_DEG)
        vz0 = REENTRY_SPEED * math.sin(gamma)
        T_ana = (vz0 + math.sqrt(vz0**2 + 2 * g0 * INITIAL_ALTITUDE)) / g0

        print(f"\nAnalytic time of flight: {T_ana:.4f} s")
        print(f"Numeric  time of flight: {T_num:.4f} s")
        print(f"Landing time error     : {abs(T_num - T_ana):.4f} s")

        # === Plot: Altitude vs Downrange ===
        plt.figure(figsize=(10, 6))
        plt.plot(x_ana, z_ana, 'k--', label='Analytic (vacuum)')
        plt.plot(x_num, z_num, 'b-',  label='Numeric (RK4, Cd=0)')
        plt.xlabel("Downrange Distance [m]")
        plt.ylabel("Altitude [m]")
        plt.title("Vacuum Trajectory: Analytic vs. RK4")
        plt.grid(True)
        plt.legend()
        plt.show()

        # === Plot 2: Altitude vs Time ===
        plt.figure(figsize=(10, 6))
        plt.plot(time_arr, z_ana, label='Analytic Altitude')
        plt.plot(time_arr, z_num, label='Numeric Altitude')
        plt.xlabel("Time [s]")
        plt.ylabel("Altitude [m]")
        plt.title("Altitude vs Time (Vacuum Comparison)")
        plt.grid(True)
        plt.legend()
        plt.show()


    # Interpolate if the last altitude is below zero
    if state_arr[-1][2] < 0:
        t0, state0, idx = interpolate_to_zero_altitude(time_arr, state_arr)
        # Replace last entry with interpolated zero-altitude state
        time_arr = np.concatenate([time_arr[:idx+1], [t0]])
        state_arr = np.concatenate([state_arr[:idx+1], [state0]])

    # Plotting results
    x_arr = state_arr[:, 0]
    y_arr = state_arr[:, 1]
    z_arr = state_arr[:, 2]
    vx_arr = state_arr[:, 3]
    vy_arr = state_arr[:, 4]
    vz_arr = state_arr[:, 5]
    speed_arr = np.sqrt(vx_arr**2 + vy_arr**2 + vz_arr**2)
    
    plt.figure(figsize=(10, 6))
    plt.plot(x_arr, z_arr, 'b-', label='Trajectory')
    plt.xlabel("Downrange Distance [m]")
    plt.ylabel("Altitude [m]")
    plt.title("Re-entry Trajectory (3DOF)")
    plt.grid(True)
    plt.legend()
    
    plt.figure(figsize=(10, 6))
    plt.plot(time_arr, speed_arr, 'r-', label='Speed')
    plt.xlabel("Time [s]")
    plt.ylabel("Speed [m/s]")
    plt.title("Vehicle Speed vs. Time")
    plt.grid(True)
    plt.legend()
    
    from mpl_toolkits.mplot3d import Axes3D
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(x_arr, y_arr, z_arr, 'g-', label='3D Trajectory')
    ax.set_xlabel("Downrange Distance [m]")
    ax.set_ylabel("Crossrange Distance [m]")
    ax.set_zlabel("Altitude [m]")
    ax.set_title("3D Trajectory of Re-entry Vehicle")
    plt.legend()
    plt.show()

    
    # Save trajectory data to a CSV file and print to console
    trajectory_data = pd.DataFrame({
        'Time [s]': time_arr,
        'x [m]': x_arr,
        'y [m]': y_arr,
        'Altitude [m]': z_arr,
        'vx [m/s]': vx_arr,
        'vy [m/s]': vy_arr,
        'vz [m/s]': vz_arr,
        'Speed [m/s]': speed_arr
    })
    trajectory_data.to_csv(OUTPUT_TRAJECTORY_CSV, index=False)
    print("Trajectory data saved to:", OUTPUT_TRAJECTORY_CSV)
    print(trajectory_data)

if __name__ == '__main__':
    # Use argparse to allow command-line options (like enabling debug prints)
    parser = argparse.ArgumentParser(description="Re-entry Trajectory Simulation")
    parser.add_argument('--debug', action='store_true', help='Enable debug prints')
    args = parser.parse_args()
    if args.debug:
        DEBUG = True
    main()
