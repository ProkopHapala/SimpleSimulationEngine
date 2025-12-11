# Tactics of Space Warfare

> This chapter is based on the original text `Latex/Taktika.tex` and provides the conceptual foundation for tactical models of how weapons, propulsion and defensive systems are used in space combat.

## 1. Nature of space engagements

The nature of combat in space follows from basic physics and from the conditions of interplanetary space. The key features are:

- **Enormous distances**
- **Low thrust of realistic engines**
- **High relative velocities**
- **(Almost) perfect visibility**

Because of these, space battles are slow to set up and to a large extent **deterministic**. Most of the "battle" actually happens in the **preparatory phase** – planning trajectories, accumulating energy, configuring the ship and the fleet – while the final, culminating engagement is very short (seconds to minutes).

Near‑perfect visibility in space usually makes it obvious which side has the upper hand. The clearly stronger side seeks contact; the weaker side tries to avoid it. If the stronger side is truly capable of catching the weaker, a battle becomes inevitable.

In a typical engagement, the ships approach each other along their trajectories until they reach the closest point of approach, and then they recede again. At closest approach, the minimum distance can range from tens of kilometres (very tight pass) up to hundreds of thousands or millions of kilometres. This point of closest approach is at the same time the point of **maximum intensity of combat and maximum effectiveness of weapons**.

Ships carry weapons usable at very different ranges. Short‑range weapons usually offer much higher firepower than long‑range ones (see also [02_Weapons](02_Weapons.md)).

The doctrine that drives the development of space armament aims to:

- **maximize firepower** in the short instant when ships are within effective range,
- **maximize effective range** of weapons.

A space battle is the culmination of weeks or months (sometimes longer) of maneuvering. Its very sharp climax often plays out in a time window of just a few seconds.

A large fraction of a warship’s equipment does not serve direct fire, but **preparation for the engagement**:

- accumulation of energy (capacitors, superconducting coils, nuclear isomers),
- preparation and reconfiguration of shields and defensive systems (see [03_Defense](03_Defense.md)),
- orientation of the hull and weapon mounts towards the enemy (internal maneuvers).

## 2. Long‑range and close‑range combat

Consider a very fast spacecraft moving at roughly $100\,\text{km/s}$. Such a ship crosses the radius of a geostationary orbit in about 6 minutes and the Earth–Moon distance in roughly one hour.

Some weapons (e.g. X‑ray lasers) can engage targets quite precisely at distances $>10^5\,\text{km}$. Other weapons (kinetic projectiles) are limited to much shorter ranges (typically $<10^3\,\text{km}$). At the same time, the destructive effect of available laser weapons is many times lower than that of heavy kinetic weapons with large mass and momentum of the projectile.

During approach, combat gradually shifts between two extremes:

- **Long‑range combat**  
  During the long approach phase, which may last tens of minutes, hours or days (for example in a pursuit), ships exchange fire using long‑range weapons. These systems are volumetrically expensive (XFEL needs a linac with undulators, other lasers require large apertures, high‑velocity kinetic guns require long acceleration tracks) and can only be mounted on the largest ships. Their effect against well‑protected targets is relatively modest – light projectiles and lasers struggle to reliably penetrate shields. Since mutual positions and orientations of ships change only slowly in this phase, shields can be configured to protect them effectively (see [03_Defense](03_Defense.md)).  
  **Effective hits** in this regime are relatively rare and the doctrine favours **large ships** with powerful long‑range weapons, strong shields and massive power sources.

- **Close‑range combat**  
  Once the distance shrinks (typically to $10^2–10^4\,\text{km}$, exceptionally less), ships can employ much more effective and cheaper weapons: primarily heavier kinetic projectiles, in some cases fired with higher cadence from many shorter accelerators. Heavier projectiles punch more easily through thin shields; high rate of fire raises the probability of hitting critical subsystems and quickly erodes protective layers. In the closest phase of the battle, **broadside salvos** become important: large numbers of projectiles are fired sideways at speeds much lower than the mutual closing speed of the ships. Heavy nuclear warheads are often delivered in this way.

Close‑range combat usually takes place in a **very short time interval**:

- for two ships passing head‑on: $\sim 0.1–10\,\text{s}$,
- in a tail chase: up to several minutes.

These times are generally shorter than what is needed for substantial **reconfiguration and reorientation of the ship** (internal maneuvers). A ship therefore has very limited ability to adapt to the situation during the close‑range phase itself, and it is crucial that it **enters the close‑range zone already in an optimal configuration**.

## 3. Maneuvering

Interplanetary spacecraft gain their speed long before the engagement itself. High‑specific‑impulse engines used to traverse interplanetary distances typically provide very low thrust ($\ll 1\,g$), and reaching cruise speeds of order $100\,\text{km/s}$ is a matter of days to weeks. The same holds for braking.

Under these conditions, the main task of piloting is to:

- put the ship or fleet on a **suitable trajectory** that either maximally approaches the target, or conversely maximally avoids the pursuing enemy,
- optimize the use of **available $\Delta v$** with regard to possible future engagements.

Key parameters are:

- $\Delta v$ – total maneuvering capability of the ship (see [04_Propulsion](04_Propulsion.md)),
- $\bar{a}$ – long‑term average acceleration (important for interception and escape),
- $a_{\text{max}}$ – peak acceleration (important for short evasive or surprise maneuvers).

In broad terms:

- a ship with **higher long‑term acceleration** $\bar{a}$ can, over time,
  - catch an enemy,
  - or avoid it if it is the one being pursued;
- **peak acceleration** $a_{\text{max}}$ determines how much a ship can "surprise" the opponent with a sudden maneuver (e.g. at the end of a chase or just before entering the close‑range zone).

## 4. Bluff and information warfare

In space, the **number** and approximate **size** of targets are usually easy to detect, so one can estimate the **overall size and rough strength** of the enemy forces reasonably well. In contrast, estimating the **exact performance and capabilities** of individual ships is difficult – each major warship is unique and its detailed technical parameters are tightly guarded.

Secrecy about true parameters (especially

- maximum thrust and $a_{\text{max}}$,
- propellant reserves and achievable $\Delta v$,
- exact performance of weapons and shields)

has enormous importance in combat. It is one of the few ways to **disrupt otherwise largely deterministic planning algorithms** on both sides.

A typical example is hiding the true maximum thrust and the duration for which it can be sustained:

- a ship may pretend to have weaker engines for a long time and only just before the engagement briefly switch to an overthrust mode,
- this may allow it to exit the enemy’s effective range if the pursuing ship cannot generate comparable thrust.

Such a maneuver is, however, **very expensive in $(\Delta v)$** and usually leaves the ship with severely limited options for the rest of the campaign. The decision whether to spend this "trick" depends critically on the technical parameters of both ships. Keeping those parameters secret greatly complicates the enemy’s planning.

## 5. Internal and external maneuvers

For ground warfare we distinguish between maneuvers that move the centre of mass of a unit (movement of a tank, marching infantry) and maneuvers that mainly change formation or orientation without significantly changing position (formation changes, rotating a turret, raising a shield).

For space combat we introduce a similar distinction:

- **External maneuver** – a maneuver that moves the centre of mass of a ship or the entire formation (changing orbit, speed, flight direction); it consumes $(\Delta v)$ and propellant.
- **Internal maneuver** – a maneuver that changes the shape, orientation or configuration of the ship without significantly moving its centre of mass (rotating segments, moving shields, reorienting weapons, redistributing masses).

External maneuvers are limited by engine thrust and propellant supply (see [04_Propulsion](04_Propulsion.md)) and consume a large amount of $(\Delta v)$. Commanders therefore try to minimize them and carry out as many operations as possible using **internal maneuvers**.

Internal maneuvers are used mainly for:

- **pointing weapons** towards the target,
- **orienting shields** so that they protect critical parts of the ship in the expected directions of incoming fire,
- **thrust vectoring** – changing the direction of the net thrust force with respect to the centre of mass.

Typical means include **rotating rings** or massive ship segments mounted in magnetic bearings, which act as counterweights when reorienting other parts. Due to conservation of angular momentum, any rotation of one part of the ship must be compensated by an opposite rotation of some other part.

In very long linear structures (e.g. particle or kinetic accelerators) fine pointing corrections can be achieved by **flexible segments** forming a kind of "spine" of the ship that allows slight curvature of the entire structure.

## 6. Ship reflexes

At the moment of closest approach and highest intensity of combat the situation is no longer manageable by a human pilot or fleet admiral. Even the main onboard computer may be too slow and too complex for some tasks.

A significant share of tasks is therefore taken over by specialized **fast feedback loops** ("ship reflexes"):

- evasive internal maneuvers,
- targeting and fire‑control corrections,
- active defense systems (laser and particle counter‑systems, interceptor projectiles – see [02_Weapons](02_Weapons.md) and [03_Defense](03_Defense.md)).

These loops are relatively simple, local and optimized for **minimal reaction time**. Above them sit more complex algorithms which, on a longer time scale (seconds to minutes), adjust their parameters, target priorities and assumptions about how the battle is likely to unfold.

A large part of battle management consists of repeated **simulations and predictions in quasi‑real time**, similar to a chess program:

- simulations use deterministic models of motion and weapon effects,
- and are continuously corrected for unexpected changes (damage, hits, unexpected enemy maneuvers).

Tuning these algorithms, choosing suitable parameters, priorities and assumptions about how the encounter will develop is a non‑trivial problem that must be solved **ad hoc** for each engagement. It requires high intelligence and expertise from the crew, who act as "configurators of reflexes" before the clash itself – during the climax, human intervention plays almost no role.

## 7. Broadside fire and inertial projectiles

Unlike on Earth, where the problem is to achieve high velocities, in space the limiting factor is more often **acceleration**. Relative velocities of ships that have been accelerating for a long time with high‑specific‑impulse engines can be comparable to or even higher than practical muzzle velocities of projectiles accelerated over limited barrel lengths.

The drive to maximize the range and effectiveness of kinetic projectiles leads to a need for **very long acceleration paths** (hundreds of metres to kilometres). Even large ships therefore usually carry only one or a few main kinetic weapons with limited rate of fire, which reduces total firepower and hit probability.

A solution is to use the **mutual relative velocity of the fleets** itself as the collision velocity of projectiles and targets:

- ships do **not fire projectiles forward** at full muzzle velocity, but release them sideways at relatively small lateral speed,
- the shot mainly serves to correct the projectile’s path relative to the tangent of the parent ship’s trajectory,
- the resulting relative projectile–target velocity is then dominated by the closing speed of the ships themselves.

With well‑designed pursuing trajectories, the required **lateral launch velocities** can be quite small. In that case it is possible to:

- release a **large number of projectiles in a short time** into a narrow bundle of trajectories,
- create a "cloud" of projectiles which:
  - ensures a high probability of hit,
  - saturates the enemy’s defensive systems,
  - allows multiple hits on the same target.

Some warships specialized for this style of combat are largely made up of **stores of projectiles** that they can rapidly release. They behave like a **directed and carefully coordinated shrapnel explosion**. A portion of the projectiles may also carry **nuclear warheads** – thanks to the absence of onboard propulsion and low launch accelerations, the construction of such warheads is simpler and cheaper (see [02_Weapons](02_Weapons.md)).

Nuclear warheads are an effective means to **cripple and strip away shields** over a large area ("soft‑killers"), so that kinetic projectiles have a clear path to the main body of the ship ("hard‑killers"). The nuclear explosions themselves usually have too little net momentum and insufficiently concentrated energy to efficiently penetrate deep into hardened targets – their main role is to disrupt defensive systems and cause widespread surface damage.

*Note:* This chapter serves as a qualitative foundation. In the file `07_Physical_model_and_equations.md` each section will be complemented by concrete mathematical models (e.g. geometry of approach, flight times, hit probabilities, optimization of broadside salvo trajectories) and by potential figures and tables stored in the `images/` directory.
