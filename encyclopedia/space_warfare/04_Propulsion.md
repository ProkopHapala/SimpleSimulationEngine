# Propulsion of spacecraft

> This chapter is based on `Latex/Pohon.tex`. It describes general requirements for military space propulsion, the division into self‑impulsive and externally‑impulsive systems, and gives an overview of the main engine types. Tactical context: [01_Tactics](01_Tactics.md). Link to ship design: [05_Ship_construction](05_Ship_construction.md). Physical models: [07_Physical_model_and_equations](07_Physical_model_and_equations.md).

## 1. Introduction and typical flight times

Propulsion of spacecraft, especially warships, must meet mutually conflicting requirements:

- high **total \(\Delta v\)** for interplanetary transfers,
- high **instantaneous acceleration** \(a\) for maneuvering in combat,
- low propellant mass fraction relative to ship mass,
- low mass of the engine itself,
- the ability to **throttle thrust and consumption** over a wide range,
- minimal demand for high‑grade energy (e.g. electricity),
- minimal waste‑heat load on the ship.

The original text refers to two key illustrations (future figures in `images/`):

- **Typical flight times in the Solar System** (`TimeOfFlight.pdf`) – the dependence of travel time on achievable velocity.
- **Specific power vs. required acceleration** (`PowerLimit.pdf`) – how required power per unit ship mass grows with demanded \(a\) and exhaust velocity.

These diagrams will later be redrawn into `images/` and analysed quantitatively in `07_Physical_model_and_equations.md`.

## 2. General propulsion options

Space engines are, by the very nature of the environment, **notoriously bad at satisfying all requirements at once**. Therefore:

- there is no single "universal" engine,
- different roles (main warships, cargo transport, missiles, infrastructure) require **different technologies**.

In this setting, the main line of "mainstream" military propulsion is **nuclear fusion with a magnetic plasma exhaust**. The era of purely chemical and purely fission thermal engines and various solar sails is treated as a transitional period before the mastery of stable fusion.

In other areas (civil transport, mining, missiles), other technologies dominate.

## 3. Self‑impulsive vs. externally‑impulsive propulsion

### 3.1 Self‑impulsive propulsion

Self‑impulsive propulsion generates impulse by **expelling some of the ship’s own mass** – classic rocketry.

- advantage: **universality** (works in vacuum everywhere),
- drawback: **must carry and consume propellant**.

Conservation of momentum says that momentum is split between the ship and the reaction mass. Because kinetic energy scales with the square of velocity, it is energetically most efficient if the medium you are pushing against is **very massive** compared to the ship (externally‑impulsive propulsion). For rockets, the opposite is true – the reaction mass is part of the ship.

Consequences:

- we want **maximum specific impulse** (as high exhaust velocity as practical),
- but the higher the **exhaust velocity**, the larger the fraction of total energy ends up in the exhaust instead of in the ship’s kinetic energy.

Example from the original text:

- we want to accelerate a ship to \(300\,\text{km/s}\) using an ion engine with \(v_e \approx 0.1c = 30\,000\,\text{km/s}\),
- the energy given to the propellant must be \(\sim 100\times\) higher than the ideal minimum (the case of an infinitely massive "reaction body").

### 3.2 Externally‑impulsive propulsion

Externally‑impulsive propulsion transfers momentum **via the surrounding medium** or other bodies:

- interaction with planetary and solar magnetic fields,
- interaction with the solar wind (plasma),
- interaction with photons (solar sails),
- gravitational maneuvers (slingshots).

Advantages:

- in principle **"infinite \(\Delta v\)"** – not limited by carried propellant,
- very energy‑efficient (we "bounce off" effectively infinite‑mass environments).

Drawbacks:

- dependence on the **external environment** (must be within a magnetosphere, plasma, near a body, etc.),
- technological complexity,
- often limited maximum force and therefore acceleration.

More detailed discussion of externally‑impulsive systems follows in chapter 6.

---

## 4. Self‑impulsive propulsion – thermal thrusters

### 4.1 Conventional thermal thrusters

Classic thermal rockets, conceptually similar to early spaceflight rockets, still have a unique ratio of **instantaneous thrust to engine mass**:

- mechanically simple, lightweight devices,
- for a short time they can expel a large mass flow of reaction material.

Drawbacks:

- limited specific impulse (set by maximum working temperature and propellant properties),
- to achieve high \(\Delta v\) they require **very large propellant reserves**.

Therefore they are mainly suitable for:

- specific tasks where short‑term high thrust is crucial,
- or where "local" material can be used as propellant.

The original text gives examples:

- **Moving asteroids and comets** – using the body’s own material as propellant for long‑term gentle maneuvers (more civil/cargo than combat application).
- **Short‑term ship maneuvers** – chemical thrusters as emergency or one‑off "trump cards" for sudden velocity changes.
- **Short‑range missiles** – brief, intense maneuvers to outsmart the parent ship or target (see [02_Weapons](02_Weapons.md)).

### 4.2 Heat sources for thermal thrusters

The original LaTeX text distinguishes several heat sources:

- **Chemical thermal thrusters**  
  - heat from chemical reactions (e.g. beryllium + hydrogen‑rich polymers),  
  - limited energy density of the fuel,  
  - advantage: very high thrust over short durations.

- **Nuclear thermal thrusters**  
  - heat source: fission reactions in fuel (uranium, plutonium, possibly californium),  
  - working temperatures of several thousand °C are possible,  
  - limit: melting point and chemical resistance of nozzle and reactor materials (tungsten, carbon/graphite have the highest melting points but are still limited).

- **Solar / electric thermal thrusters**  
  - solar concentrators or electric heaters,  
  - more suited for less urgent, energy‑efficient maneuvers.

Nuclear thermal thrusters also raise the issue of **neutron radiation** and shielding of the crew compartment, which is a major structural challenge (see [05_Ship_construction](05_Ship_construction.md)).

---

## 5. Self‑impulsive propulsion – advanced concepts

The original chapter lists several advanced concepts, only sketched for now:

### 5.1 Nuclear pulse propulsion

- analogous to Project Orion – the ship is accelerated by a series of small nuclear explosions,
- advantages: potentially very high \(\Delta v\) and thrust,
- drawbacks: enormous structural loads, radiation exposure, political and environmental issues,
- more suitable for large ships or heavy tugs than for finely maneuvering warships.

### 5.2 Continuous fusion magnetic nozzles

The text only names several variants:

- **magnetic mirror lines**,  
- **toroidal magnetic confinement** (tokamak / stellarator),  
- **Polywell‑type devices**.

General idea:

- fusion plasma in a magnetic confinement device serves as both energy source and working medium,
- plasma is exhausted through a magnetic nozzle to produce thrust,
- limits of plasma confinement, density and temperature (tokamak physics) become central.

### 5.3 Thermal plasma thrusters and EM acceleration

- **Thermal plasma thrusters** – similar to nuclear thermal, but with plasma as the working medium.  
- **Electromagnetic acceleration of ions and plasma** – various ion engines and plasma thrusters (VASIMR‑like).  
- **Electromagnetic acceleration of solid mass** – railguns and coilguns used as ship propulsion (mass‑drivers).  
- **Direct fusion propulsion** – using fusion reaction products directly as exhaust.

In the model these systems will be handled mainly through:

- specific impulse,  
- maximum achievable thrust,  
- \(P/m\) ratio (power per unit mass of the propulsion system),  
- thermal losses and radiator requirements.

---

## 6. Externally‑impulsive propulsion

### 6.1 Solar sails

- use pressure of solar radiation on a large reflective sail,  
- advantages: no carried propellant, in principle inexhaustible \(\Delta v\),  
- drawbacks: very low acceleration, dependence on distance from the star.

Military use:

- mainly for long‑term transfers and "stealthy" motion of small objects,  
- only limited use in ship‑to‑ship combat where high maneuverability is needed.

### 6.2 Magnetic and electrostatic sails

- **Magnetic sails** – large current loops generating a magnetic field that interacts with the solar wind,  
- **Electrostatic sails** – large charged structures interacting with ions in the plasma.

In both cases:

- the propulsion system is "spread out" and has a very large physical size,  
- more suitable for lighter craft and infrastructure than for heavy warships,  
- in warfare, more of a **strategic tool** (e.g. gradual repositioning or towing of large objects).

### 6.3 Electrodynamic tether systems

Electrodynamic tether systems (conductors in a magnetosphere driven by current) can:

- generate forces by interacting with the planetary magnetic field,
- be used for braking or accelerating in orbit,
- form "tracks" for space transport (see also infrastructure in `Pohon.tex` and the section on space transport infrastructure).

Militarily:

- a key part of space infrastructure (stations, tugs, tracks),
- warships can dock to them, but in "deep" space they must rely on their own engines.

---

## 7. Space transport infrastructure and war

For civil and commercial transport, the following play major roles:

- **fixed magnetic tracks and loops** around major bodies,
- **braking and acceleration tracks** tied to orbital stations,
- **space elevators and orbital tugs**.

From a military point of view:

- well‑developed infrastructure is one of the main factors of **economic and strategic power**,  
- and represents **key targets** for attack and defence (see [06_Scenarios_and_history](06_Scenarios_and_history.md) for possible narrative scenes).

---

## 8. Link to simulations

When building `07_Physical_model_and_equations.md`, this chapter will provide in particular:

- relations between \(I_{sp}\), \(\Delta v\), propellant mass and ship mass,
- the relation between propulsion system power \(P\), exhaust velocity and achievable thrust,
- models of long‑term average acceleration \(\bar{a}\) vs. short‑term \(a_{\text{max}}\),
- constraints imposed by thermal management and radiator area,
- typical flight times between key Solar System locations for various propulsion types.

These models will serve as a basis for simulating tactical scenarios in [01_Tactics](01_Tactics.md) and for estimating maneuverability of ships in the planned game/simulation.
