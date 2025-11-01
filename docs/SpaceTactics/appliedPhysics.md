# `appliedPhysics.h`

This header file provides a collection of fundamental physical constants and inline functions for various physics calculations used throughout the simulation. It serves as a centralized library for common physics formulas.

## Constants

A wide range of physical and time-related constants are defined for use in calculations.

### Physical Constants

*   `const_GravAccel`: Standard gravity acceleration (9.81 m/s²).
*   `const_Graviational`: The gravitational constant (6.6743015e-11 N·(m/kg)²).
*   `const_Rgas`: The ideal gas constant (8.314... J/(mol·K)).
*   `const_eV`: Electron-volt in Joules.
*   `const_ElectronCharge`: The elementary charge in Coulombs.
*   `const_Bonltzman`: The Boltzmann constant (1.380649e-23 J/K).
*   `const_massProton`: The mass of a proton in kg.
*   `const_massElectron`: The mass of an electron in kg.
*   `const_AU`: The astronomical unit in meters (149.6e9 m).
*   `const_SolarRadEarth`: Solar radiation intensity at Earth's orbit (1366.1 W/m²).

### Energy Density Constants

*   `const_EkgAnihil`: Energy from matter-antimatter annihilation per kg (in Joules).
*   `const_EkgDT`: Energy from Deuterium-Tritium fusion per kg.
*   `const_EkgU`: Energy from Uranium fission per kg.
*   `const_EktTNT`: Energy equivalent of a kiloton of TNT in Joules.
*   `const_EMWy`: Energy equivalent of a Megawatt-year in Joules.

### Other Constants

*   `const_heatCapacityRatio_...`: Heat capacity ratios (γ) for monoatomic, diatomic, and water vapor gases.
*   `const_hour`, `const_day`, `const_month`, `const_year`: Time durations in seconds.

## Functions

*   `char* timeInfo(char* s, double t_sec)`: Converts a time in seconds to a human-readable string (e.g., "1.5 hours", "2.3 years").
*   `double solarRadDist_SI(double r)`: Calculates the solar radiation intensity in W/m² at a given distance `r` from the sun.
*   `double kineticEnergy(double v, double mass)`: Calculates the kinetic energy of an object.
*   `double gunEnergy(...)`: A model for calculating the muzzle energy of a propellant-based gun, considering factors like powder mass and chamber volume.
*   `double armorThicnessFactor_MomentumModel(...)`: A model to calculate the effective thickness of sloped armor based on momentum refraction upon impact.
*   `Vec3d centralGravityForce(const Vec3d& d, double Mm)`: Calculates the gravitational force vector between two bodies, given the displacement vector `d` and the product of their masses `Mm`.
*   `double tsielkovsky_speed(double payload, double vexh)`: Implements the Tsiolkovsky rocket equation to calculate the maximum change in velocity (delta-v) based on payload fraction and exhaust velocity.
*   `double tsielkovsky_payload(double deltaV, double vexh)`: The inverse of the rocket equation, calculating the required payload fraction to achieve a given `deltaV`.
*   `double jetEfficiency(double expansionRatio, double kappa)`: Calculates the thermodynamic efficiency of a rocket nozzle based on the gas expansion ratio.
*   `double exhaustVelocity(double T, double molarMass, ...)`: Calculates the exhaust velocity of a rocket engine based on temperature, molar mass of the propellant, and engine efficiency.
