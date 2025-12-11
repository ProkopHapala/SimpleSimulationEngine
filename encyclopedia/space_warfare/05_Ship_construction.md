# Construction of spacecraft

> This chapter is based on `Latex/konstrukce_lode.tex` and on notes in `Texts/Design_kosmickych_lodi.md`. It describes the basic layout of warships, the advantages of large ships, the difference between linear and circular designs, and the main structural challenges. For propulsion physics see [04_Propulsion](04_Propulsion.md), for tactics see [01_Tactics](01_Tactics.md).

## 1. Why large ships?

As discussed in other chapters, there are many reasons why **large spacecraft are preferred** over greater numbers of small ships:

1. **Scaling laws of engines**  
   Engines of larger physical size (especially pulsed and fusion nuclear engines) are generally more efficient – they offer a better power‑to‑mass ratio.

2. **Large weapons**  
   Effective long‑range weapons – linear and circular accelerators, free‑electron lasers, large apertures – require long paths and large supporting structures.

3. **Internal maneuvers**  
   Physically connected modules (by tethers or structural members) allow configuration changes by internal maneuvers (without shifting the centre of mass) instead of expensive external maneuvers (see [01_Tactics](01_Tactics.md)).

4. **Shared and redundant systems**  
   An interconnected network of sub‑units can share highly redundant systems (reactors, radiators, sensors, shields). The overall resilience of the "super‑ship" is higher than that of the same number of isolated small vessels.

These advantages are counterbalanced by serious **engineering difficulties**, described below.

## 2. Main structural challenges

### 2.1 Sparsely filled volume

- It is undesirable for vital, vulnerable systems to occupy a compact, densely filled volume – that would increase the probability of a crippling hit.
- Many subsystems themselves require sparse volumes:
  - electromagnetic tether systems,
  - solar sails,
  - radiator panels,
  - low‑pressure plasma vessels of fusion reactors,
  - long accelerator lines.

The result is a **large but sparse structure** – something between a bridge, a crane and a stretched spider web.

### 2.2 Structural integrity

It must be ensured that individual parts of the ship **do not move unpredictably** relative to each other even under forces and accelerations during maneuvers:

- long trusses and tethers must carry tension and bending in low gravity but with significant inertial forces,
- linear dimensions may far exceed those of terrestrial cranes or bridges.

This puts strong demands on:

- lightweight yet strong **composite materials**,  
- **tensegrity structures** (combinations of struts in compression and cables in tension),  
- active damping and control of vibrations.

### 2.3 "Locomotor system" for internal maneuvers

Maximizing the rate at which the ship can reconfigure its shape requires a kind of **"locomotor system"**:

- an analogue of bones and muscles in living organisms,
- actuators (electric motors, linear drives) coupled to trusses and tethers,
- use of **rotating rings** and other massive parts as counterweights.

The goals are:

- to use existing functional components (reactors, tanks, weapons) as "muscles" and "masses",  
- to minimize extra dead mass dedicated solely to maneuvering,  
- to avoid over‑stressing the structure with impulsive forces.

### 2.4 Armor and protection

Distributing systems over a large volume makes **effective armoring** more difficult:

- the simplest solution would be to armor all critical systems under a single protective shell,
- but in practice they are scattered along the length of the ship and in different directions.

In addition:

- protective shields are usually **highly directional** (they cover only a narrow solid angle),  
- they must be coordinated with internal maneuvers and ship rotation (see [03_Defense](03_Defense.md)).

## 3. Principle of "more tethers, fewer beams"

One of the basic design principles is to **use cables instead of solid beams wherever possible**:

- a taut cable is lighter than a solid beam carrying the same tensile load,
- beams are used mainly in compression and for local stiffening.

Another principle is to use **centrifugal force** to tension cables:

- rotating rings and frames create tension without the need for heavy rigid members,
- this is particularly suitable for circular ship designs.

This reliance on flexible, low‑stiffness structures and on inertial forces associated with rotation leads to very complex **dynamics of ship motion**, which must be:

- carefully simulated during design,  
- continuously evaluated by onboard computers during maneuvers.

## 4. Linear ship concepts

### 4.1 Basic scheme

Linear ship concepts start from the fact that:

- main weapons (linear accelerators, free‑electron lasers, railguns) are **straight and long**,
- main thrust forces act along the **longitudinal axis** of the ship,
- a linear shape minimizes frontal cross‑section in head‑on combat.

Advantages:

- minimal frontal area when the main weapons are pointed directly at a distant enemy ⇒ minimum exposed target area,
- a well‑defined **"forward combat direction"**, around which shields and weapon mounts can be optimised,
- simpler geometry for some weapon types (e.g. linear accelerators).

Drawbacks:

- less opportunity to tension the structure centrifugally,
- less natural arrangement for some externally‑impulsive methods (large‑radius magnetic loops).

### 4.2 Tactical implications

A linear ship is an ideal platform for:

- **front‑on long‑range engagements** (small frontal profile, strong shield "cap" in a single direction),
- concentration of fire into a narrow cone (see [01_Tactics](01_Tactics.md), long‑range combat).

Conversely, in close‑range environments or when attacked from multiple directions it may struggle to:

- rapidly redirect its main weapons to another azimuth,
- reconfigure shields to cover multiple directions at once.

## 5. Circular ship concepts

### 5.1 Basic scheme

Circular ship concepts have complementary pros and cons compared to linear ones:

- the primary structural element is a **ring** (or several rings),
- the ring can simultaneously serve as:
  - an **acceleration track** (circular accelerator for projectiles or particles),
  - a **load‑bearing structure** tensioned by centrifugal force,
  - potentially as a structural component of a magnetic sail.

Advantages:

- for a given radius the ship can be **lighter** (the structure works more in tension than in bending),
- the **circular track** can be used to store energy and projectiles,
- projectiles can be launched over a wide range of azimuths along the circumference.

Drawbacks:

- to achieve the same projectile acceleration as a linear barrel, a **longer path** is required (roughly $\pi$ times longer than a straight section),
- in typical long‑range combat only a small solid angle is tactically relevant – the ability to "fire anywhere" is not fully used.

### 5.2 Tactical implications

A circular ship can excel in:

- **close passes and fleet fly‑bys**, where rapid fire into many azimuths without large internal maneuvers is valuable,
- situations requiring **very fast vectoring of fire** all around the ship.

On the other hand:

- in purely head‑on long‑range duels it may not offer an advantage over a linear ship,
- integration of some propulsion types is more complex (see [04_Propulsion](04_Propulsion.md)) if the main engine must be aligned with the centre of mass.

### 5.3 Relation to magnetic sails

The original text notes that:

- an ideal **magnetic sail** typically spans hundreds to thousands of kilometres,
- a ship ring of order ~10 km is much smaller and far too heavy,
- as a magnetic sail it is therefore unsuitable – small rings would require unreasonably large currents and fields.

Magnetic sails thus often form **a separate extended system**, only weakly coupled to the ship (by tethers, fields) and not identical with the main load‑bearing structure.

## 6. Lightweight structures of large size (outline)

The original LaTeX chapter ends here with only a heading; for further development it is useful to plan for:

- **truss and tensegrity structures** – three‑dimensional frameworks and tensegrity designs,
- **modular segments** – the ability to reconfigure the ship by docking and undocking modules,
- **active damping** – sensor networks measuring vibrations and actuator networks damping them.

This section will later be expanded and linked to specific models in `07_Physical_model_and_equations.md` (stresses in cable structures during maneuvers, limits on maximum accelerations, resonance frequencies).

## 7. Locomotor system for zero‑$\Delta v$ maneuvers (outline)

The second unfinished subsection of the original text will cover in particular:

- **mass distribution** (reactors, tanks, weapons, radiators) as potential counterweights,
- a **network of actuators** capable of transferring forces and torques between segments,
- **acceleration limits** for sensitive components (crew, electronics) during fast internal maneuvers.

Together with [01_Tactics](01_Tactics.md) and [03_Defense](03_Defense.md) this should lead to a model that tells us:

- how fast a ship can realistically change its orientation and configuration,
- how this depends on its shape (linear vs. circular) and mass distribution,
- what are the "safe" and "risky" regimes of internal maneuvering.

## 8. Links to other chapters and simulations

- **Propulsion** – location of engines, thrust vectoring, link between structure and maximum achievable acceleration: [04_Propulsion](04_Propulsion.md).
- **Tactics** – use of ship shape and internal maneuvering in combat: [01_Tactics](01_Tactics.md).
- **Defense** – ship shape and distribution of shields/armor: [03_Defense](03_Defense.md).

In `07_Physical_model_and_equations.md` this chapter will feed into:

- models of stresses and allowed accelerations for long structures,
- simple models of dynamics of rotating rings and tethers,
- estimates of structural mass fractions (structure mass vs. payload) for different ship concepts.
