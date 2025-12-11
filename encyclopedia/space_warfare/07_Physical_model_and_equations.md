# Physical model and equations

> This chapter collects physical models, simplified equations and parameters needed for simulating space combat and designing spacecraft. It is based on [01_Tactics](01_Tactics.md), [02_Weapons](02_Weapons.md), [03_Defense](03_Defense.md), [04_Propulsion](04_Propulsion.md), [05_Ship_construction](05_Ship_construction.md) and notes in the `Texts/` directory.

The goal is to have **one central place** where:

- all physical quantities and parameters are defined,
- the used approximations and simplifications are recorded,
- code for a simulation or game can be implemented directly from here.

## 1. Orbital mechanics and ship kinematics

### 1.1 Basic quantities

- $m$ – ship mass [kg]
- $v$ – ship velocity [m/s]
- $\Delta v$ – total "maneuver budget" in velocity [m/s]
- $a(t)$ – instantaneous acceleration [m/s²]
- $\bar{a}$ – long‑term average acceleration [m/s²]
- $a_{\text{max}}$ – short‑term peak acceleration [m/s²]

### 1.2 Two‑body approximation

For a first, coarse model of motion between bodies in the Solar System we can use a simplified two‑body approximation (ship + central body). In this document we focus more on **time scales and $\Delta v$) requirements** than on precise orbital mechanics.

For combat simulations (time scales minutes to days) it is often sufficient to:

- model ships in an **inertial frame** with given initial velocities,
- add small $\Delta v$ kicks from thrusters.

Detailed Keplerian dynamics can be left to external libraries; here it is enough to know typical orders of magnitude:

- cruise velocities ($\sim 10^2\,\text{km/s}$),
- campaign $\Delta v$ ($\sim 10–10^3\,\text{km/s}$).

### 1.3 Time to contact and closest approach (sketch)

For two objects with initial relative distance $L_0$ and relative velocity $v_{\text{rel}}$ (neglecting gravity):

$$
 t_{\text{contact}} \approx \frac{L_0}{v_{\text{rel}}}.
$$

For combat (see [01_Tactics](01_Tactics.md)) we care mainly about:

- the window when the target is within **effective weapon range**,
- how fast the target can change direction ($a_{\text{max}}$) to avoid hits.

A more detailed model (including acceleration during approach) can be added later.

## 2. Propulsion and energetics

Based on [04_Propulsion](04_Propulsion.md) and related notes.

### 2.1 Rocket equation (self‑impulsive propulsion)

For engines that consume their own reaction mass, the Tsiolkovsky rocket equation applies:

$$
 \Delta v = v_e \ln \left(\frac{m_0}{m_f}\right),
$$

where:

- $v_e$ – effective exhaust velocity [m/s],
- $m_0$ – initial mass (ship + propellant + payload) [kg],
- $m_f$ – final mass (ship + payload) [kg].

Specific impulse $I_{sp}$ is defined as:

$$
 I_{sp} = \frac{v_e}{g_0},
$$

with $g_0 \approx 9.81\,\text{m/s}^2$.

### 2.2 Power, thrust and exhaust velocity

For an engine with exhaust velocity $v_e$, mass flow $\dot{m}$ and thrust $T$:

$$
 T = \dot{m} \, v_e.
$$

Mechanical power in the jet:

$$
 P_{\text{jet}} = \frac{1}{2} \dot{m} v_e^2 = \frac{T^2}{2 \dot{m}}.
$$

If electrical / nuclear power $P$ is available with efficiency $\eta$, then roughly

$$
 P_{\text{jet}} \leq \eta P.
$$

This leads to a trade‑off:

- higher $v_e$ ⇒ lower propellant consumption, but for given $P$ smaller thrust,  
- lower $v_e$ ⇒ higher thrust, but lower $I_{sp}$ and poorer $\Delta v$.

### 2.3 Long‑term vs. short‑term acceleration

We distinguish:

- $\bar{a}$ – acceleration that the ship can sustain for long periods (limited by power and $\Delta v$),
- $a_{\text{max}}$ – short‑term peak acceleration (limited by structure and short boost capability).

In simulations we will often use a simplified model where each ship is described by $\Delta v$, $\bar{a}$ and $a_{\text{max}}$, from which chase / escape options are derived (see [01_Tactics](01_Tactics.md)).

## 3. Propagation of laser beams

Based on [02_Weapons](02_Weapons.md) and `Vykony_kosmickych_lodi.md`.

### 3.1 Diffraction limit

For a circular aperture of diameter $D$ and wavelength $\lambda$, the minimum angular width of the main lobe is approximately

$$
 \theta \approx 1.22\, \frac{\lambda}{D}.
$$

Spot diameter at range $L$:

$$
 d \approx \theta L \approx 1.22 \, \frac{\lambda L}{D}.
$$

In practical units (approximation from the notes):

$$
 \lambda_{[\mu m]} \cdot L_{[\text{thousand km}]} \approx D_{[\text{m}]} \cdot d_{[\text{m}]}.
$$

### 3.2 Intensity at target

Total laser power $P_L$, neglecting losses, is distributed over area $A = \pi (d/2)^2$. Average intensity:

$$
 I \approx \frac{P_L}{A} = \frac{4 P_L}{\pi d^2}.
$$

To damage material, $I$ must exceed some threshold (depending on exposure time, material properties and cooling). This can be modelled either by

- a simple threshold $I > I_{\text{crit}}$, or
- a more sophisticated **erosion model** (see chapter 6 below).

## 4. Particle beams

### 4.1 Divergence of charged beams (sketch)

For strong charged‑particle beams, collective effects are complex. For our purposes we introduce a **smooth** parametrised attenuation rather than a sudden cut‑off at some range.

We define $L_{\text{eff}}$ as a characteristic length scale over which the current density significantly decreases, e.g. via exponential or power‑law attenuation. A simple on‑axis model:

$$
 j(L) = j_0 \exp\left[-\left(\frac{L}{L_{\text{eff}}}\right)^\alpha\right],
$$

where $\alpha \sim 1{-}2$ is a shape parameter.

The parameter $L_{\text{eff}}$ lumps together the effects of

- particle type (electron, proton, heavy ion),
- average energy and current density,
- focusing quality (optics, magnetic lenses),
- **Coulomb self‑repulsion** in the beam,
- **collisional / diffusive scattering** in the environment (plasma, residual atmosphere).

Future refinements may explicitly distinguish between

- purely Coulomb‑driven beam expansion (faster, accelerated spread),
- diffusive scattering (slower growth, $\propto \sqrt{t}$).

#### Coulombic expansion – simple sketch

Consider a beam of particles (charge $q$, mass $m$) with current $I$ and transverse radius $r_b$. The space‑charge field can be very roughly represented by an effective transverse acceleration

$$
 a_\perp \sim \frac{q E_\perp}{m} \propto \frac{I}{m v_z r_b},
$$

where $v_z$ is the longitudinal velocity. The narrower the beam (smaller $r_b$) and the higher the current $I$, the larger the repulsive force and the faster the beam "bloats".

A simplified equation for the growth of beam radius along the axis (with $z = v_z t$) is

$$
 \frac{d^2 r_b}{dz^2} \approx \frac{1}{v_z^2} a_\perp(r_b) \propto \frac{I}{m v_z^3 r_b}.
$$

This leads to accelerated growth $r_b(z)$ (faster than $\sqrt{z}$), whereas a purely diffusive process would give $r_b \propto \sqrt{t} \propto \sqrt{z}$. For game‑level modelling, detailed coefficients can be hidden in $L_{\text{eff}}$, but the direction of dependencies (higher current $\Rightarrow$ faster loss of focus, more relativistic beams $\Rightarrow$ slower expansion) should agree with these relations.

### 4.2 Penetration and damage

For a simple model we distinguish between

- **surface damage** (small penetration depth – typical for sub‑relativistic charged beams),
- **volume damage** (larger penetration depth – neutrons, very relativistic particles).

Stopping power and range can in a first approximation be represented by parameters

- $R$ – characteristic penetration depth for a given weapon / material,  
- $k$ – coefficient converting deposited energy into "damage".

More detailed models can later be added using simple tables.

## 5. Ballistics of kinetic projectiles

Based on [02_Weapons](02_Weapons.md) and `Vykony_kosmickych_lodi.md`.

### 5.1 Acceleration by force

For constant driving force $F$ and projectile mass $m$:

$$
 a = \frac{F}{m}, \quad v(x) = \sqrt{2 a x} = \sqrt{\frac{2 F x}{m}}.
$$

Exit velocity $v$ grows with the **square root of path length** and **square root of $F/m$**.

### 5.2 Power‑limited acceleration (railguns)

From the notes:

- time of flight of a projectile over distance $D$:

$$
 T(v) = \frac{D}{v} + T_{\text{accel}}(v).
$$

If input power $P$ is limited and projectile kinetic energy is $E = \tfrac{1}{2} m v^2$, then

$$
 T_{\text{accel}} \approx \frac{E}{P} = \frac{m v^2}{2P}.
$$

Total time:

$$
 T(v) = \frac{D}{v} + \frac{m v^2}{2P}.
$$

Minimising $T(v)$ gives an optimal

$$
 v_{\text{opt}} \propto \left(\frac{P D}{m}\right)^{1/3}.
$$

At this optimum, **acceleration time and cruise time are comparable**.

### 5.3 Target deviation under maneuvering

If the target maneuvers during the projectile’s flight with typical acceleration $g$ (e.g. $g \sim a_{\text{max}}$), it will deviate from the predicted path by

$$
 d_s = \frac{1}{2} g T^2,
$$

where $T$ is projectile flight time (including acceleration). In practice, the diameter of the "lethal cloud" of projectiles (e.g. a broadside salvo) must exceed $d_s$, otherwise the target can evade.

## 6. Shield and defence model

Based on [03_Defense](03_Defense.md).

### 6.1 Erosion of thin shields

A simple model for **thin shields** (foils, plasma curtains):

- define surface mass density $\sigma$ [kg/m²],
- each hit with energy $E$ over area $A$ reduces $\sigma$ by $\Delta \sigma \propto E/A$,
- once $\sigma < \sigma_{\text{min}}$, the shield is perforated at that spot.

### 6.2 Probability of stopping projectiles

For active defence (lasers, particle beams, interceptor projectiles) we can define:

- $R_{\text{det}}$ – maximum detection range,
- $R_{\text{eng}}$ – effective engagement range of the defence system,
- $f$ – firing rate [shots/s],
- $p_{\text{hit}}$ – probability that a single attempt destroys / deflects one projectile.

For an incoming salvo of projectiles the penetration probability can be approximated as

$$
 P_{\text{penetrate}} \approx \exp(- N_{\text{eff}} p_{\text{hit}}),
$$

where $N_{\text{eff}}$ is the effective number of engagement attempts (depends on the time the projectiles stay in range and on $f$).

## 7. Nuclear weapons (0D model)

This section summarizes a very simplified (0‑dimensional) model of fission and fusion warheads that allows us to estimate **burn fraction** and **yield** without solving full neutron transport.

### 7.1 Fission bomb – neutron population growth

Basic parameters:

- $M_f$ – mass of fissile material (e.g. $^{235}\text{U}$, $^{239}\text{Pu}$) [kg],
- $\rho$ – density of the active core [kg/m³],
- $k_{\text{eff}}$ – effective multiplication factor (after all leakage and losses),
- $\Lambda$ – generation time for neutrons [s],
- $t_{\text{dis}}$ – time until the core expands mechanically so much that the chain reaction stops [s].

The number of neutrons (and approximately the number of fissions in each generation) behaves as a geometric series:

$$
 N_{n+1} = k_{\text{eff}} N_n, \quad N_n = N_0 k_{\text{eff}}^{n},
$$

with generation index $n$. During time $t_{\text{dis}}$ roughly

$$
 n_{\text{gen}} \approx \frac{t_{\text{dis}}}{\Lambda}
$$
generations occur.

The total number of fissions during the explosion is approximated by

$$
 N_{\text{fiss}} \approx N_0 \frac{k_{\text{eff}}^{n_{\text{gen}}} - 1}{k_{\text{eff}} - 1},
$$

capped above by the total number of available fissile nuclei $N_{\text{tot}}$. The **burn fraction** is then

$$
 f_{\text{burn}} = \min\left(1, \, \frac{N_{\text{fiss}}}{N_{\text{tot}}}\right).
$$

The yield (neglecting losses into radiation leakage, case kinetic energy, etc.) is

$$
 Y \approx f_{\text{burn}} \, N_{\text{tot}} \, E_{\text{fiss}},
$$

where $E_{\text{fiss}}$ is energy released per fission ($\sim 200\,\text{MeV}$). In practice, $f_{\text{burn}}$ is often tens of percent, not 100\,\%.

The critical mass $M_{\text{crit}}$ for a given material and geometry depends on $\rho$ and cross‑sections; roughly, increasing density reduces $M_{\text{crit}}$ approximately as

$$
 M_{\text{crit}} \propto \frac{1}{\rho^2},
$$

reflecting the fact that compressing the core (e.g. in implosion designs) greatly reduces the needed mass.

### 7.2 Thermonuclear warheads and Lawson criterion

For the thermonuclear stage (e.g. D–T mixture) the key quantity is the product of density, temperature and confinement time. The **Lawson criterion** is

$$
 n T \tau_E \ge (n T \tau)_\text{crit},
$$

where

- $n$ – particle density of the fuel [m⁻³],
- $T$ – plasma temperature [keV],
- $\tau_E$ – energy confinement time [s],
- $(n T \tau)_\text{crit}$ – critical value depending on the reaction (for D–T of order $10^{21}$ $[{keV·s·m}^{-3}]$ ).

In a 0D model we assume that during confinement time $\tau$ a certain **burn fraction** $f_{\text{burn,DT}}$ is achieved, which can be approximated as a rising function of

$$
 \chi = \frac{n T \tau}{(n T \tau)_\text{crit}}.
$$

A simple saturating form is, for example,

$$
 f_{\text{burn,DT}} \approx 1 - \exp(-C_\chi \, \chi),
$$

with $C_\chi \sim 1$ encoding detailed losses. The thermonuclear yield is then

$$
 Y_{\text{fus}} \approx f_{\text{burn,DT}} \, N_{\text{DT}} \, E_{\text{fus}},
$$

where $N_{\text{DT}}$ is the number of D–T pairs and $E_{\text{fus}}$ the energy per fusion reaction.

Total yield is the sum of fission and fusion parts:

$$
 Y_{\text{tot}} \approx Y + Y_{\text{fus}}.
$$

These 0D models allow us to vary parameters (fuel mass, compression, implosion quality) in a simulation and see their effect on yield without solving full transport.

## 8. Structural mechanics of ships (outline)

Based on [05_Ship_construction](05_Ship_construction.md).

Basic estimates:

- tension in a tether for acceleration $a$: $T = m_{\text{attached}} a$,
- maximum allowed acceleration for tether material with tensile strength $\sigma_{\text{max}}$, cross‑section $S$ and attached mass $m$:

$$
 a_{\text{max}} \approx \frac{\sigma_{\text{max}} S}{m}.
$$

Similarly, we can estimate stresses in long trusses during rotation and internal maneuvers.

More detailed models (including resonance frequencies and damping) can be added later if the simulation needs to limit internal maneuver rates.

## 9. Summary and use in simulation

This file is intended to serve as a **working whitepaper** for implementing:

- strategic / tactical simulations of space combat,
- simple calculators of ship and weapon performance.

Next steps:

- add concrete numerical examples and plots (stored in `images/`),
- decide what level of realism is appropriate for the game/simulation,
- for each equation, specify the **parameter range** where it makes sense and where it breaks down.
