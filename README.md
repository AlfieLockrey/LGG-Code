# LGG Code

# Inputs for Riad 12.7mm 8g Charge Comparison
<div>
# COMBUSTION CHAMBER VALUES
D_c = 30e-3                     # Diameter of the Combustion Chamber
L0_c = 15e-3                    # Length of the Combustion Chamber
P0_c = 100e3                    # Initial pressure in the Combustion Chamber (before detonation)
C = 8e-3                       # Charge Mass 5 to 20 g

gamma_c = 1.2238                # γ Gamma for the combustion products

BR_exp = 0.81837                # BURN RATE Exponent
u1 = 3.102e-8                   # Burning Rate Constant
e1 = 1.27e-4                    # Propellant half web size - assumes hollow cylindrical powder
R1 = u1 / e1                    # Experimentally determined Burn Rate Coefficient
density_propel = 1.6e3          # δ Propellant Density
ForceConst_propel = 1.159e6     # λ Propellant Force Constant
CoVolume_propel = 0.8e-3        # η Propellant Co-Volum
V0_c = L0_c * math.pi * (D_c / 2)**2

# PISTON VALUES
D_pis = 30e-3            # Diameter of piston
m_pis = 110e-3          # Mass of the piston
mu_static_pis = 0       # Coefficient of friction for the piston against the pump tube
mu_dynamic_pis = 0
allowPistonRearward = False

# PUMP TUBE VALUES
P0_pt = 10e6       # Initial pressure in the pump tube (ahead of the piston)
L0_pt = 1.3        # Length of the pump tube
D_pt = 30e-3        # Diameter of the pump tube
gamma_lg = 1.66     # Gamma for the light gas

# RUPTURE DISK VALUES
P_rupt = 35e+6  # Pressure at which the rupture disk ruptures

# BARREL VALUES
L_b = 1.9           # Length of the barrel
D_b = 12.7e-3       # Diameter of the barrel
P0_b = 133          # Initial upstream pressure in the barrel

# IMPACT CHAMBER
gamma_ic = 1.4      # Gamma for the (near vaccum) gas in the impact chamber
V_ic = 1

# PROJECTILE & SABOT VALUES
m_pr = 17e-3    # Projectile mass
m_sb = 0        # Sabot Mass0.47e-3
mu_sb = 0     # Sabot Friction coeff.

# SIMULATION DETAILS
delta_t = 1e-6  # Time step length

</div>