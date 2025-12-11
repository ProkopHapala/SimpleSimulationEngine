# Defense of Spacecraft

> This chapter is based on the original text `Latex/Obrana.tex`. It introduces the basic classes of defensive measures and their links to the other chapters: [01_Tactics](01_Tactics.md), [02_Weapons](02_Weapons.md) and [04_Propulsion](04_Propulsion.md).

## 1. Basic types of defense

In this framework, the defense of spacecraft can be divided into several basic mechanisms:

1. **Evasion without shifting the centre of mass**  
   Fast internal maneuvers that change the orientation and configuration of the ship, but not its overall trajectory.

2. **Evasion with a shift of the centre of mass**  
   External maneuvers using the main engines – changing the trajectory of the ship (consuming $ \Delta v$).

3. **Thin shields**  
   Lightweight, replaceable protective layers (foils, plasma or aerosol clouds, dust layers), effective mainly against beam weapons and small kinetic projectiles.

4. **Massive shields**  
   Thick layers of material or entire space objects (asteroids, moons, planets) that absorb or disperse extreme amounts of energy.

In the following text this classification is refined and linked to specific types of attack from [02_Weapons](02_Weapons.md).

## 2. Principles of passive protection

### 2.1 Hiding behind massive obstacles

One of the most effective defensive methods is to **hide behind a sufficiently massive obstacle** such that penetrating or destroying it would require the attacker to expend an extreme amount of energy and resources – more than is economical compared to the target itself.

Typical obstacles:

- natural bodies: planets, moons, asteroids,
- large artificial structures: massive stations, stockpiles of raw material.

In a realistic setting **there is no such thing as a "sufficiently strong material" or a universal "energy shield"** that would reliably withstand modern weapons, especially nuclear and high‑energy kinetic weapons. Passive protection always works in terms of **probabilities of damage** and economic balance (how much mass and energy the attacker must invest).

### 2.2 Thin shields and deflectors

Against lasers and charged‑particle beams, **thin shields** are advantageous:

- thin metallic or composite foils,
- plasma and aerosol clouds,
- systems of "deflectors" and mirrors that reflect or scatter incoming beams.

Typically they:

- protect only a **limited solid angle** around the ship,
- when correctly oriented, allow **very high local effectiveness** at modest mass cost,
- can be rapidly deployed and retracted (internal maneuvers – see [01_Tactics](01_Tactics.md)).

## 3. Principles of active protection

The second major group of defensive techniques is **active deflection or destruction of incoming threats**:

- deflecting beams or particle streams (with magnetic fields, plasma),
- deflecting or breaking up projectiles (lasers, particle beams, counter‑projectiles),
- jamming the guidance of missiles.

This group is broad and each principle is suited to a different type of offensive weapon. A basic classification (for details see [02_Weapons](02_Weapons.md)):

- **Against energy beams**  
  deflectors, plasma lenses, aerosol or dust clouds, specially tuned optical materials.

- **Against charged‑particle beams**  
  magnetic and electric fields, plasma envelopes that scatter and deflect the beam.

- **Against neutral beams and hard radiation**  
  options are more limited; one often relies on massive shielding material and increasing the effective size of the target (spreading the deposited energy over a large volume).

- **Against kinetic projectiles**  
  laser and particle systems that ablate surface layers and deflect the trajectory, or fragment the projectile; also "cloud armor" – clouds of particles or screening structures.

- **Against guided missiles**  
  combinations of sensor jamming, laser and particle fire, and anti‑missiles.

## 4. Evasion – internal and external maneuvers

The original text distinguishes two types of "evasion":

1. **Without shifting the centre of mass** – high‑quality internal maneuvers:
   - rapid rotation of the hull so that the smallest possible profile is presented towards the direction of attack,
   - moving thin shields into a narrow solid angle from which fire is expected,
   - rearranging masses (propellant tanks, cargo) so that they shield sensitive systems.

2. **With a shift of the centre of mass** – a true evasive maneuver:
   - changing the trajectory using main engines,
   - paid for in $\Delta v$ and limitations on options in later phases of the campaign (see [04_Propulsion](04_Propulsion.md)).

As a rule, it is desirable to **minimize the number of external maneuvers**, which are expensive in propellant. Most defensive actions should be handled by proper ship design (layout of mass, shields and weapons) and by internal maneuvers (see also [01_Tactics](01_Tactics.md), the section on internal maneuvers).

## 5. Realistic limits of "shields"

From the viewpoint of physics and materials it is important to emphasize:

- there is no universal "force field" resistant to all weapon types,
- each shield is **optimized for a certain class of threats** (thin foils against lasers, plasma clouds against particle beams, massive armor against kinetics),
- practical defense is always a **combination** of several principles:
  - geometry (hiding behind bodies),
  - material shields,
  - active systems (lasers, particle beams, interceptor projectiles),
  - ship maneuvers.

In practice, ships are designed so that they can:

- **concentrate protection** into a narrow cone from which the most intense fire is expected,
- and at the same time reconfigure this protection on time scales of seconds to minutes as the situation evolves.

## 6. Links to other chapters

- The tactical context of defense (when to hide behind a body, when to actively cover incoming fire, when to spend $\Delta v$ to break contact) is discussed in detail in [01_Tactics](01_Tactics.md).
- The nature of individual threats and physical parameters of weapons are described in [02_Weapons](02_Weapons.md).
- Structural and material aspects of shields and supporting structures will be covered in [05_Ship_construction](05_Ship_construction.md).
- Quantitative models (probability of shield penetration, erosion of thin layers, required material thickness for a given type of hit) will be collected in [07_Physical_model_and_equations](07_Physical_model_and_equations.md).
