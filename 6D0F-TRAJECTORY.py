import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Function to call the Runge-Kutta Numerical Integration Factors into the Script
def rkfactors(nrkstages):
    # Function to call the Runge-Kutta Numerical Integration Factors into the Script
    # Enforce integer and positive input
    if not isinstance(nrkstages, int) or nrkstages <= 0:
        raise ValueError("nrkstages must be a positive integer")
    if nrkstages > 5:
        raise ValueError("Please enter an integer value of between 1 and 5 for the number of RungeKutta Numerical Integration Stages. It is recommended that 4 stages be used. ")
    frk = [0]*nrkstages
    if nrkstages == 1:
        frk[0] = 1.0
    elif nrkstages == 2:
        frk[0] = 0.5
        frk[1] = 1.0
    elif nrkstages == 3:
        frk[0] = 0.6
        frk[1] = 0.6
        frk[2] = 1.0
    elif nrkstages == 4:
        frk[0] = 0.25
        frk[1] = 1/3
        frk[2] = 0.5
        frk[3] = 1.0
    elif nrkstages == 5:
        frk[0] = 0.25
        frk[1] = 1/6
        frk[2] = 3/8
        frk[3] = 0.5
        frk[4] = 1.0
    return frk

# Function to determine free stream atmospheric properties, dependent on altitude
def ATMOSPHERICPROPERTIESFUNCTION(Z):
    IDEAL_GAS_CONSTANT_AIR = 0.287053
    IDEAL_GAS_SPECIFIC_HEAT_RATIO_AIR = 1.4
    if Z <= 11000:
        TEMPERATURE = 288.15 - 0.0065 * Z
        PRESSURE = 101.325 * ((TEMPERATURE / 288.15) ** 5.256)
    elif (Z > 11000) and (Z <= 25000):
        TEMPERATURE = 216.65
        PRESSURE = 22.65 * math.exp(1.73 - 0.000157 * Z)
    else:
        TEMPERATURE = 141.94 + 0.00299 * Z
        PRESSURE = 2.488 * ((TEMPERATURE / 216.65) ** -11.388)
    RHO = PRESSURE / (IDEAL_GAS_CONSTANT_AIR * TEMPERATURE)
    SOUND = math.sqrt(IDEAL_GAS_SPECIFIC_HEAT_RATIO_AIR * 1000 * IDEAL_GAS_CONSTANT_AIR * TEMPERATURE)
    return RHO, SOUND

# Function to determine the gravitational force components in the body-fixed reference frame
def GRAVITYFUNCTION(xg, yg, zg, psi, theta, phi):
    GRAVITY_SEA_LEVEL = 9.80665
    MEAN_EARTH_RADIUS = 63566766
    GRAVITY_X = (-GRAVITY_SEA_LEVEL / MEAN_EARTH_RADIUS) * (xg * math.cos(psi) * math.cos(theta) + yg * math.sin(psi) * math.cos(theta) - (MEAN_EARTH_RADIUS - 2 * zg) * math.sin(theta))
    GRAVITY_Y = (-GRAVITY_SEA_LEVEL / MEAN_EARTH_RADIUS) * (xg * (math.cos(psi) * math.sin(theta) * math.sin(phi) - math.sin(psi) * math.cos(phi)) + yg * (math.sin(psi) * math.sin(theta) * math.sin(phi) + math.cos(psi) * math.cos(phi)) + (MEAN_EARTH_RADIUS - 2 * zg) * math.cos(theta) * math.sin(phi))
    GRAVITY_Z = (-GRAVITY_SEA_LEVEL / MEAN_EARTH_RADIUS) * (xg * (math.cos(psi) * math.sin(theta) * math.cos(phi) + math.sin(psi) * math.sin(phi)) + yg * (math.sin(psi) * math.sin(theta) * math.cos(phi) - math.cos(psi) * math.sin(phi)) + (MEAN_EARTH_RADIUS - 2 * zg) * math.cos(theta) * math.cos(phi))
    return GRAVITY_X, GRAVITY_Y, GRAVITY_Z

# Projectile Data and Initial Conditions - User Inputs
# 6-DOF Trajectory Prediction Code – User Defined Projectile

# User Inputs - Projectile Data
MASS = 43.7  # [kg]
REFERENCE_LENGTH = 0.155  # [m]
REFERENCE_AREA = 0.25 * math.pi * REFERENCE_LENGTH**2  # [m^2]
IXX = 0.1444  # [kg.m^2]
IYY = 1.7323  # [kg.m^2]
IZZ = 1.7323  # [kg.m^2]

# Call AeroCoefficients Data
# Using pandas to read variables from the file.
# Assumes the file is in CSV format with columns matching the variable names.
data = pd.read_csv('DENEL 155mm ASSEGAI M2000 SHELL AEROCOEFFICIENTS DATA.csv')
MACH_data = data['MACH_data'].to_numpy()
CDRAGF_LINEAR_data = data['CDRAGF_LINEAR_data'].to_numpy()
CDRAGF_QUADRATIC_data = data['CDRAGF_QUADRATIC_data'].to_numpy()
CLIFTF_LINEAR_data = data['CLIFTF_LINEAR_data'].to_numpy()
CLIFTF_CUBIC_data = data['CLIFTF_CUBIC_data'].to_numpy()
CMAGNUSF_COEFF_data = data['CMAGNUSF_COEFF_data'].to_numpy()
CDAMPINGF_COEFF_data = data['CDAMPINGF_COEFF_data'].to_numpy()
CSPINM_COEFF_data = data['CSPINM_COEFF_data'].to_numpy()
COVERTURNINGM_COEFF_data = data['COVERTURNINGM_COEFF_data'].to_numpy()
CMAGNUSM_COEFF_data = data['CMAGNUSM_COEFF_data'].to_numpy()
CDAMPINGM_COEFF_data = data['CDAMPINGM_COEFF_data'].to_numpy()

# User Inputs - Initial Conditions
LAUNCH_VELOCITY = 481  # [m/s]
LAUNCH_ANGLE = 67.5  # [degrees]
INITIAL_HEIGHT = 0  # [m]
INITIAL_SPIN_RATE = 847.2  # [rad/s]

# User Inputs - Numerical Method
TIMESTEP = 0.0001  # [s]
NRKSTAGES = 4

# Calculations, Assigning Array Variables and Initial Conditions – Automatic Code
# Calculations - Projectile Data and Initial Conditions
LAUNCH_ANGLE_RADS = LAUNCH_ANGLE * math.pi / 180
INITIAL_X_VELOCITY = LAUNCH_VELOCITY * math.cos(LAUNCH_ANGLE_RADS)
INITIAL_Z_VELOCITY = LAUNCH_VELOCITY * math.sin(LAUNCH_ANGLE_RADS)
INIT_ANGULAR_SPIN = INITIAL_SPIN_RATE
INIT_YAW_RADS = -LAUNCH_ANGLE_RADS

# Initialise Arrays (using lists)
time = [0]
x = [0]
y = [0]
z = [INITIAL_HEIGHT]
vxe = [INITIAL_X_VELOCITY]
vye = [0]
vze = [INITIAL_Z_VELOCITY]
ve = [LAUNCH_VELOCITY]
phi = [0]
theta = [INIT_YAW_RADS]
psi = [0]
vx = []
vy = []
vz = []
v = []
Mach = []
omegax = []
omegay = []
omegaz = []
omega = []
yaw = []
yawdot = []

# Assign initial conditions to arrays
# Compute initial body-axis velocities
vx_val = vxe[0] * math.cos(psi[0]) * math.cos(theta[0]) + vye[0] * math.sin(psi[0]) * math.cos(theta[0]) - vze[0] * math.sin(theta[0])
vy_val = (vxe[0] * (math.cos(psi[0]) * math.sin(theta[0]) * math.sin(phi[0]) - math.sin(psi[0]) * math.cos(phi[0])) +
          vye[0] * (math.sin(psi[0]) * math.sin(theta[0]) * math.sin(phi[0]) + math.cos(psi[0]) * math.cos(phi[0])) +
          vze[0] * math.cos(theta[0]) * math.sin(phi[0]))
vz_val = (vxe[0] * (math.cos(psi[0]) * math.sin(theta[0]) * math.cos(phi[0]) + math.sin(psi[0]) * math.sin(phi[0])) +
          vye[0] * (math.sin(psi[0]) * math.sin(theta[0]) * math.cos(phi[0]) - math.cos(psi[0]) * math.sin(phi[0])) +
          vze[0] * math.cos(theta[0]) * math.cos(phi[0]))
vx.append(vx_val)
vy.append(vy_val)
vz.append(vz_val)
v.append(math.sqrt(vx_val**2 + vy_val**2 + vz_val**2))
_, SOUND = ATMOSPHERICPROPERTIESFUNCTION(z[0])
Mach.append(ve[0] / SOUND)
omegax.append(INIT_ANGULAR_SPIN)
omegay.append(0)
omegaz.append(0)
omega.append(math.sqrt(omegax[0]**2 + omegay[0]**2 + omegaz[0]**2))
yaw.append(math.acos(vx[0] / v[0]))
yawdot.append(0)

# Call Runge-Kutta factors
frk = rkfactors(NRKSTAGES)

# Perform Runge-Kutta Numerical Integration of Differential Equations of Motion
it = 0  # Set it variable to first iteration (Python indexing starts at 0)
while (z[it] > 0) or (time[it] == 0):  # Loop runs until projectile impacts the ground
    it += 1
    U = vx[it - 1]
    V = vy[it - 1]
    W = vz[it - 1]
    TRANSVEL = math.sqrt(U**2 + V**2 + W**2)
    MACH_current = Mach[it - 1]
    P = omegax[it - 1]
    Q = omegay[it - 1]
    R = omegaz[it - 1]
    ROTVEL = math.sqrt(P**2 + Q**2 + R**2)
    X = x[it - 1]
    Y = y[it - 1]
    Z = z[it - 1]
    THETAX = phi[it - 1]
    THETAY = theta[it - 1]
    THETAZ = psi[it - 1]
    ATTACK = yaw[it - 1]
    ATTACKDOT = yawdot[it - 1]
    
    RHO_current, _ = ATMOSPHERICPROPERTIESFUNCTION(Z)
    AERO_FORCE_CONSTANT = 0.5 * RHO_current * REFERENCE_AREA / MASS
    AERO_MOMENT_CONSTANT = 0.5 * RHO_current * REFERENCE_AREA * REFERENCE_LENGTH
    
    CDRAGF_LINEAR = np.interp(MACH_current, MACH_data, CDRAGF_LINEAR_data)
    CDRAGF_QUADRATIC = np.interp(MACH_current, MACH_data, CDRAGF_QUADRATIC_data)
    CLIFTF_LINEAR = np.interp(MACH_current, MACH_data, CLIFTF_LINEAR_data)
    CLIFTF_CUBIC = np.interp(MACH_current, MACH_data, CLIFTF_CUBIC_data)
    CMAGNUSF_COEFF = np.interp(MACH_current, MACH_data, CMAGNUSF_COEFF_data)
    CDAMPINGF_COEFF = np.interp(MACH_current, MACH_data, CDAMPINGF_COEFF_data)
    COVERTURNINGM_COEFF = np.interp(MACH_current, MACH_data, COVERTURNINGM_COEFF_data)
    CMAGNUSM_COEFF = np.interp(MACH_current, MACH_data, CMAGNUSM_COEFF_data)
    CSPINM_COEFF = np.interp(MACH_current, MACH_data, CSPINM_COEFF_data)
    CDAMPINGM_COEFF = np.interp(MACH_current, MACH_data, CDAMPINGM_COEFF_data)
    
    CDRAGF = CDRAGF_LINEAR + CDRAGF_QUADRATIC * (math.sin(ATTACK) ** 2)
    CLIFTF = CLIFTF_LINEAR * math.sin(ATTACK) + CLIFTF_CUBIC * (math.sin(ATTACK) ** 3)
    CMAGNUSF = CMAGNUSF_COEFF * P * math.sin(ATTACK) * REFERENCE_LENGTH / TRANSVEL
    CDAMPINGF = CDAMPINGF_COEFF * ATTACKDOT * REFERENCE_LENGTH / TRANSVEL
    COVERTURNINGM = COVERTURNINGM_COEFF * math.sin(ATTACK)
    CMAGNUSM = CMAGNUSM_COEFF * P * math.sin(ATTACK) * REFERENCE_LENGTH / TRANSVEL
    CSPINM = CSPINM_COEFF * P * REFERENCE_LENGTH / TRANSVEL
    CDAMPINGM = CDAMPINGM_COEFF * ATTACKDOT * REFERENCE_LENGTH / TRANSVEL
    
    for RKSTAGE in range(1, NRKSTAGES + 1):
        GRAVITY_X, GRAVITY_Y, GRAVITY_Z = GRAVITYFUNCTION(X, Y, Z, THETAZ, THETAY, THETAX)
        THRUSTF = THRUSTFORCEFUNCTION(time[it - 1])
        if ROTVEL == 0:
            vx_new = vx[it - 1] + frk[RKSTAGE - 1] * TIMESTEP * (GRAVITY_X + THRUSTF / MASS + AERO_FORCE_CONSTANT * (CLIFTF * (V ** 2 + W ** 2) - CDRAGF * TRANSVEL * U))
            vy_new = vy[it - 1] + frk[RKSTAGE - 1] * TIMESTEP * (GRAVITY_Y + AERO_FORCE_CONSTANT * (CMAGNUSF * TRANSVEL * W - CDRAGF * TRANSVEL * V - CLIFTF * U * V))
            vz_new = vz[it - 1] + frk[RKSTAGE - 1] * TIMESTEP * (GRAVITY_Z + AERO_FORCE_CONSTANT * (- CDRAGF * TRANSVEL * W - CMAGNUSF * TRANSVEL * V - CLIFTF * U * W))
            v_new = math.sqrt(vx_new ** 2 + vy_new ** 2 + vz_new ** 2)
            omega_x_new = omegax[it - 1] + frk[RKSTAGE - 1] * TIMESTEP * (AERO_MOMENT_CONSTANT * (TRANSVEL ** 2) * (1 / IXX) * CSPINM)
            omega_y_new = omegay[it - 1] + frk[RKSTAGE - 1] * TIMESTEP * (AERO_MOMENT_CONSTANT * TRANSVEL * (1 / IYY) * (COVERTURNINGM * W + CMAGNUSM * V))
            omega_z_new = omegaz[it - 1] + frk[RKSTAGE - 1] * TIMESTEP * (AERO_MOMENT_CONSTANT * TRANSVEL * (1 / IZZ) * (CMAGNUSM * W - COVERTURNINGM * V))
            omega_new = math.sqrt(omega_x_new ** 2 + omega_y_new ** 2 + omega_z_new ** 2)
        else:
            vx_new = vx[it - 1] + frk[RKSTAGE - 1] * TIMESTEP * (R * V - Q * W + GRAVITY_X + THRUSTF / MASS + AERO_FORCE_CONSTANT * (CLIFTF * (V ** 2 + W ** 2) - CDRAGF * TRANSVEL * U))
            vy_new = vy[it - 1] + frk[RKSTAGE - 1] * TIMESTEP * (P * W - R * U + GRAVITY_Y + AERO_FORCE_CONSTANT * TRANSVEL * (CMAGNUSF * W - CDRAGF * V - CLIFTF * U * V / TRANSVEL - CDAMPINGF * TRANSVEL * R / ROTVEL))
            vz_new = vz[it - 1] + frk[RKSTAGE - 1] * TIMESTEP * (Q * U - P * V + GRAVITY_Z + AERO_FORCE_CONSTANT * TRANSVEL * (CDAMPINGF * TRANSVEL * Q / ROTVEL - CDRAGF * W - CMAGNUSF * V - CLIFTF * U * W / TRANSVEL))
            v_new = math.sqrt(vx_new ** 2 + vy_new ** 2 + vz_new ** 2)
            omega_x_new = omegax[it - 1] + frk[RKSTAGE - 1] * TIMESTEP * (AERO_MOMENT_CONSTANT * (TRANSVEL ** 2) * (1 / IXX) * CSPINM - R * Q * (IZZ - IYY) / IXX)
            omega_y_new = omegay[it - 1] + frk[RKSTAGE - 1] * TIMESTEP * (AERO_MOMENT_CONSTANT * TRANSVEL * (1 / IYY) * (COVERTURNINGM * W + CMAGNUSM * V + CDAMPINGM * TRANSVEL * Q / ROTVEL) - R * P * (IXX - IZZ) / IYY)
            omega_z_new = omegaz[it - 1] + frk[RKSTAGE - 1] * TIMESTEP * (AERO_MOMENT_CONSTANT * TRANSVEL * (1 / IZZ) * (CMAGNUSM * W - COVERTURNINGM * V + CDAMPINGM * TRANSVEL * R / ROTVEL) - P * Q * (IYY - IXX) / IZZ)
            omega_new = math.sqrt(omega_x_new ** 2 + omega_y_new ** 2 + omega_z_new ** 2)
        x_new = x[it - 1] + frk[RKSTAGE - 1] * TIMESTEP * (U * math.cos(THETAY) * math.cos(THETAZ) + 
                V * (math.sin(THETAX) * math.sin(THETAY) * math.cos(THETAZ) - math.cos(THETAX) * math.sin(THETAZ)) +
                W * (math.cos(THETAX) * math.sin(THETAY) * math.cos(THETAZ) + math.sin(THETAX) * math.sin(THETAZ)))
        y_new = y[it - 1] + frk[RKSTAGE - 1] * TIMESTEP * (U * math.cos(THETAY) * math.sin(THETAZ) +
                V * (math.sin(THETAX) * math.sin(THETAY) * math.sin(THETAZ) + math.cos(THETAX) * math.cos(THETAZ)) +
                W * (math.cos(THETAX) * math.sin(THETAY) * math.sin(THETAZ) - math.sin(THETAX) * math.cos(THETAZ)))
        z_new = z[it - 1] + frk[RKSTAGE - 1] * TIMESTEP * (- U * math.sin(THETAY) + V * math.sin(THETAX) * math.cos(THETAY) + W * math.cos(THETAX) * math.cos(THETAY))
        phi_new = phi[it - 1] + frk[RKSTAGE - 1] * TIMESTEP * (P + Q * math.sin(THETAX) * math.tan(THETAY) + R * math.cos(THETAX) * math.tan(THETAY))
        theta_new = theta[it - 1] + frk[RKSTAGE - 1] * TIMESTEP * (Q * math.cos(THETAX) - R * math.sin(THETAX))
        psi_new = psi[it - 1] + frk[RKSTAGE - 1] * TIMESTEP * ((1 / math.cos(THETAY)) * (R * math.cos(THETAX) + Q * math.sin(THETAX)))
        
        # Update variables for next RK stage
        vx[it - 1] = vx_new
        vy[it - 1] = vy_new
        vz[it - 1] = vz_new
        v[it - 1] = v_new
        omegax[it - 1] = omega_x_new
        omegay[it - 1] = omega_y_new
        omegaz[it - 1] = omega_z_new
        omega[it - 1] = omega_new
        x[it - 1] = x_new
        y[it - 1] = y_new
        z[it - 1] = z_new
        phi[it - 1] = phi_new
        theta[it - 1] = theta_new
        psi[it - 1] = psi_new
        
        U = vx[it - 1]
        V = vy[it - 1]
        W = vz[it - 1]
        TRANSVEL = v[it - 1]
        P = omegax[it - 1]
        Q = omegay[it - 1]
        R = omegaz[it - 1]
        ROTVEL = omega[it - 1]
        X = x[it - 1]
        Y = y[it - 1]
        Z = z[it - 1]
        THETAX = phi[it - 1]
        THETAY = theta[it - 1]
        THETAZ = psi[it - 1]
    
    time.append(time[it - 1] + TIMESTEP)
    vxe_val = (vx[it - 1] * math.cos(theta[it - 1]) * math.cos(psi[it - 1]) +
               vy[it - 1] * (math.sin(phi[it - 1]) * math.sin(theta[it - 1]) * math.cos(psi[it - 1]) - math.cos(phi[it - 1]) * math.sin(psi[it - 1])) +
               vz[it - 1] * (math.cos(phi[it - 1]) * math.sin(theta[it - 1]) * math.cos(psi[it - 1]) + math.sin(phi[it - 1]) * math.sin(psi[it - 1])))
    vxe.append(vxe_val)
    vye_val = (vx[it - 1] * math.cos(theta[it - 1]) * math.sin(psi[it - 1]) +
               vy[it - 1] * (math.sin(phi[it - 1]) * math.sin(theta[it - 1]) * math.sin(psi[it - 1]) + math.cos(phi[it - 1]) * math.cos(psi[it - 1])) +
               vz[it - 1] * (math.cos(phi[it - 1]) * math.sin(theta[it - 1]) * math.sin(psi[it - 1]) - math.sin(phi[it - 1]) * math.cos(psi[it - 1])))
    vye.append(vye_val)
    vze_val = -vx[it - 1] * math.sin(theta[it - 1]) + vy[it - 1] * math.sin(phi[it - 1]) * math.cos(theta[it - 1]) + vz[it - 1] * math.cos(phi[it - 1]) * math.cos(theta[it - 1])
    vze.append(vze_val)
    ve_val = math.sqrt(vxe_val ** 2 + vye_val ** 2 + vze_val ** 2)
    ve.append(ve_val)
    RHO_new, SOUND_new = ATMOSPHERICPROPERTIESFUNCTION(z[it - 1])
    Mach.append(ve_val / SOUND_new)
    yaw_val = math.acos(vx[it - 1] / v[it - 1])
    yaw.append(yaw_val)
    yawdot.append(abs((yaw[-1] - yaw[-2]) / TIMESTEP) if len(yaw) > 1 else 0)
    
    vx.append(vx[it - 1])
    vy.append(vy[it - 1])
    vz.append(vz[it - 1])
    v.append(v[it - 1])
    omegax.append(omegax[it - 1])
    omegay.append(omegay[it - 1])
    omegaz.append(omegaz[it - 1])
    omega.append(omega[it - 1])
    x.append(x[it - 1])
    y.append(y[it - 1])
    z.append(z[it - 1])
    phi.append(phi[it - 1])
    theta.append(theta[it - 1])
    psi.append(psi[it - 1])

# Transpose arrays (convert to numpy arrays and then column vectors)
time = np.array(time).reshape(-1, 1)
vx = np.array(vx).reshape(-1, 1)
vy = np.array(vy).reshape(-1, 1)
vz = np.array(vz).reshape(-1, 1)
v = np.array(v).reshape(-1, 1)
Mach = np.array(Mach).reshape(-1, 1)
omegax = np.array(omegax).reshape(-1, 1)
omegay = np.array(omegay).reshape(-1, 1)
omegaz = np.array(omegaz).reshape(-1, 1)
omega = np.array(omega).reshape(-1, 1)
x = np.array(x).reshape(-1, 1)
y = np.array(y).reshape(-1, 1)
z = np.array(z).reshape(-1, 1)
vxe = np.array(vxe).reshape(-1, 1)
vye = np.array(vye).reshape(-1, 1)
vze = np.array(vze).reshape(-1, 1)
ve = np.array(ve).reshape(-1, 1)
phi = np.array(phi).reshape(-1, 1)
theta = np.array(theta).reshape(-1, 1)
psi = np.array(psi).reshape(-1, 1)
yaw = np.array(yaw).reshape(-1, 1)
yaw_degrees = yaw * 180 / math.pi
yawdot = np.array(yawdot).reshape(-1, 1)

# Trajectory Plots – Run automatically by code
plt.figure(1)
plt.plot(x, z, 'r-')
plt.xlabel(' x [m]')
plt.ylabel(' z [m]')
plt.xlim([0, x[-1, 0]])
plt.ylim([0, np.max(z)])
plt.title('Range vs Height')

plt.figure(2)
ax = plt.axes(projection='3d')
ax.plot3D(x.flatten(), y.flatten(), z.flatten(), 'k-')
ax.set_xlabel(' x [m]')
ax.set_ylabel(' y [m]')
ax.set_zlabel('z [m]')
plt.title('3D Trajectory')

plt.figure(3)
plt.plot(time, yaw_degrees, 'k-')
plt.xlabel('Time [s]')
plt.ylabel('Yaw Angle [degrees]')
plt.title('Yaw Angle vs Time')

plt.show()

# Trajectory Data Table – Run automatically by code
Trajectory_Data = pd.DataFrame({
    'time': time.flatten(),
    'x': x.flatten(),
    'y': y.flatten(),
    'z': z.flatten(),
    'vxe': vxe.flatten(),
    'vye': vye.flatten(),
    'vze': vze.flatten(),
    'v': v.flatten(),
    'Mach': Mach.flatten(),
    'omegax': omegax.flatten(),
    'yaw_degrees': yaw_degrees.flatten()
})
Trajectory_Data.Properties = {}  # Dummy attribute for compatibility with MATLAB-like syntax
Trajectory_Data = Trajectory_Data.rename(columns={
    'time': 'Time [s]',
    'x': 'Range [m]',
    'y': 'Drift [m]',
    'z': 'Height [m]',
    'vxe': 'X Velocity [m/s]',
    'vye': 'Y Velocity [m/s]',
    'vze': 'Z Velocity [m/s]',
    'v': 'Velocity [m/s]',
    'Mach': 'Mach Number',
    'omegax': 'Axial Spin Rate [m/s]',
    'yaw_degrees': 'Yaw Angle [degrees]'
})
print(Trajectory_Data)
    
