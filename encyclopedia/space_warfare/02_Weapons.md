# Weapons of Space Warfare

> This chapter is based on the original text `Latex/Zbrane.tex` and describes the main categories of weapons considered for combat in the Solar System. For tactical context see [01_Tactics](01_Tactics.md), for modelling see [07_Physical_model_and_equations](07_Physical_model_and_equations.md).

## 1. Overview of weapon categories

In this framework, space weapons are divided into:

1. **Energy beams**  
   Directed electromagnetic radiation – various kinds of lasers, and possibly high‑power EMP.
2. **Particle beams**  
   Beams of charged or neutral particles with high energy.
3. **Solid projectiles**  
   Kinetic weapons – electromagnetically accelerated projectiles, nuclear projectiles, electrically accelerated aerosol.
4. **Guided missiles**  
   Small autonomous craft with their own propulsion and usually a nuclear warhead.

In this chapter we focus on:

- physical limits of each weapon type (range, effectiveness, penetration),
- their link to tactical scenarios (long‑range vs. close‑range combat – see [01_Tactics](01_Tactics.md)),
- basic modelling ideas for future simulations.

---

## 2. Energy beams

Energy beams include lasers, directed electromagnetic pulses (EMP) and other forms of EM radiation. Common features:

- **Almost unlimited ammunition** – limited only by available energy and cooling capacity.
- **Propagation speed ≈ c** – reaction times are mostly set by pointing and electronics.

Disadvantages:

- **Extreme demands on peak power** during the shot.
- Often **large mass and size** of the device relative to delivered firepower.
- For lasers: a **diffraction limit** – the required aperture area grows with range (see also `Vykony_kosmickych_lodi.md`).

### 2.1 Gas lasers

Gas lasers use gaseous active media and are better suited than classic solid‑state lasers because of:

- high power relative to the mass of the active medium,
- the possibility of direct electrical pumping or pumping via gas dynamics.

Examples:

- **CO₂ lasers** – vibrational transitions in the IR band; less suitable due to long wavelength and large aperture requirements.
- **Excimer lasers** – shorter wavelengths, more suitable for long‑range engagement and small spot sizes.

Disadvantages of gas lasers:

- complex gas handling (cooling, circulation),
- limited firing rate due to heating of the working gas and cooling constraints,
- low efficiency (few percent), requiring very powerful electrical power sources.

### 2.2 Plasma‑dynamic and nuclear‑pumped lasers

#### 2.2.1 Plasma‑dynamic lasers

The goal is to achieve **very short wavelengths** (X‑ray lasers) using electronic transitions in highly ionized plasma of heavy elements. Typically:

- a **pinch discharge** is used to rapidly compress and heat the plasma,
- during subsequent expansion, a population inversion and laser action occur.

Problems:

- low efficiency – only a small fraction of transitions contributes to population inversion,
- need for a complex reactor or pulsed plasma source.

#### 2.2.2 Lasers pumped by nuclear explosives

The most powerful lasers can be pumped directly by **nuclear explosions**:

- a very short, extremely intense pulse,
- the device is typically **destroyed on the first shot**, which limits them to special tactical situations.

Conceptual uses:

- one‑off, extreme flashes to punch through strong defences,
- pumping large diffraction‑limited beams, ideally through a network of deflectors/satellites (see [Poznamky_k_boji_ve_vesmiru](../Texts/Poznamky_k_boji_ve_vesmiru.md)).

#### 2.2.3 Neutron‑pumped lasers

A nuclear reactor (fission or fusion) producing neutrons can secondarily pump lasers, for example as follows:

- neutrons interact with suitable nuclei (\(^3\)He, \(^6\)Li, etc.),
- high‑energy particles or non‑equilibrium electronic excitations are produced,
- these can drive laser action in a heavy gas or plasma.

This concept is technologically complex but promises a **more direct conversion of nuclear energy** into coherent radiation.

### 2.3 Free‑electron lasers (FEL)

Free‑electron lasers amplify an electromagnetic wave by decelerating relativistic electrons in a **magnetic undulator**:

- energetically expensive, but **tunable in wavelength**,
- can be driven either by a conventional accelerator or directly by a nuclear reactor (extraction of electrons or nuclear fragments).

Advantages for space warfare:

- smoothly tunable frequency – the ability to adapt to properties of enemy shields,
- potential for relatively high powers when well coupled to a reactor.

Disadvantages:

- significant mass and complexity,
- difficult integration into a compact warship (more suitable for large hub ships or fixed installations).

### 2.4 Gamma lasers and exotic concepts

Here we only sketch the possibilities of **gamma‑ray lasers based on nuclear transitions** and **electron–positron annihilation lasers**. These concepts assume:

- the existence of suitable nuclear isomers,
- the ability to safely store large amounts of annihilation fuel.

In the current framework they should be treated more as **theoretical limits** that define how "hard" radiation is even conceivable in principle.

### 2.5 Diffraction limit and beam geometry

The figure in the original text (`Difraction.pdf`) illustrates the relation between **wavelength \(\lambda\)**, **aperture \(D\)**, distance \(L\) and spot diameter \(d\) on target. Approximately:

\[
\theta \sim 1.22\, \frac{\lambda}{D}, \quad d \sim \theta L.
\]

In the notes `Vykony_kosmickych_lodi.md` this is written in the convenient form

\[
\lambda_{[\mu m]} \cdot L_{[\text{thousand km}]} \approx D_{[\text{m}]} \cdot d_{[\text{m}]}.
\]

This simple approximation will be elaborated in `07_Physical_model_and_equations.md` and illustrated with plots.

---

## 3. Particle beams

Particle beams can be divided into **charged** and **neutral** beams.

General properties:

- strong interaction with matter (ionization, nuclear reactions),
- in charged beams, strong sensitivity to electric and magnetic fields,
- in neutral beams, acceleration and focusing are technically difficult.

### 3.1 Charged‑particle beams

Types:

- electron beams,
- ion beams (light and heavy ions),
- in principle also plasma beams.

Advantages:

- relatively **simple acceleration** in electromagnetic accelerators,
- possibility of high energies (1–10 GeV, up to TeV in very large machines),
- in the target they can produce spallation reactions and secondary radioactivity.

Disadvantages:

1. **Collective spreading**  
   - like‑charged particles repel each other,  
   - high‑intensity beams tend to expand rapidly,  
   - at relativistic energies the situation improves due to time dilation and length contraction, but it remains a limiting factor.

2. **Easy defence and sensitivity to environment**  
   - beams can be bent strongly by magnetic and electric fields,  
   - they are influenced by planetary magnetospheres, the solar wind,  
   - one can build active **magnetic shields** against them.

3. **Small penetration depths** (for subrelativistic energies)  
   - charged particles interact strongly with atomic shells,  
   - they can be stopped by thin material layers (surface ablation rather than deep destruction).

#### 3.1.1 Use in active defence

Despite their limitations, **charged beams** are very attractive for **short‑range active defence**:

- at distances of units to tens of kilometres they can ablate surface layers and deflect a projectile by metres to tens of metres,
- they can be miniaturized into many barrels with high firing rate (kHz) and modest power.

Combined with powerful computers and radar/optical sensors they can form a "fire curtain" against:

- micrometeoroids,
- small kinetic projectiles,
- particles of electrically accelerated aerosol.

### 3.2 Neutral beams

Neutral beams (of neutrons or atoms) remove some fundamental drawbacks of charged beams:

- they have no net charge ⇒ much less sensitivity to external fields,
- space‑charge effects are greatly reduced ⇒ much less spontaneous spreading,
- potentially **large penetration depth** into matter.

Problems:

- **acceleration and focusing** of neutral particles is much more demanding,
- for neutrons one mostly has isotropic sources (reactors, spallation sources),
- atomic beams require first accelerating ions and then neutralizing them.

#### 3.2.1 Neutron beams

- produced in reactors or spallation sources as **hot neutrons**,  
- incoherent and nearly isotropic,  
- focusing is difficult (diffraction on crystal lattices, extreme fields),  
- cooling and re‑acceleration is very energy‑inefficient.

In the current concept there is **no practical way** known to use neutrons as a directional beam weapon, but:

- a free neutron at relativistic speed can travel Earth–Mars distances within its half‑life,
- neutrons have excellent penetration and initiate nuclear reactions deep inside the target.

#### 3.2.2 Atomic beams

Concept:

- accelerate cations of heavy elements (e.g. W, Th) together with hydrogen,  
- just before the muzzle **neutralize** the beam (electron capture),  
- in rapid adiabatic expansion the hydrogen carries away most of the heat and heavy elements condense into droplets/particles.

Challenges:

- the heat released during neutralization must not spread the beam too much,
- very fast cooling (10–100 µs) over a relatively short acceleration path (100–1000 m) is required,
- in the end, the penetration capability of atomic beams is limited in a similar way as for charged beams.

### 3.3 Electrically accelerated aerosol

A compromise between atomic beams and classical kinetic projectiles is **electrically accelerated aerosol**:

- solid dust grains of microscopic size (e.g. diamond, tungsten carbide, uranium carbide),
- the grains are charged by capturing ions/electrons and then accelerated in RF cavities,
- before leaving the muzzle they are neutralized again.

Advantages:

- simpler technology than for atomic beams,
- most of the heat released during neutralization goes into internal degrees of freedom of the particle (lattice vibrations),
- there is no critical spreading problem like in a gas beam.

Disadvantages:

- achievable exit velocities are lower than for atomic beams,
- penetration is limited – one gets micro‑craters of millimetre scale,
- flight times to target are longer than for very fast atomic or charged beams.

---

## 4. Solid projectiles

Solid projectiles are classical **kinetic weapons** – momentum and kinetic energy of the projectile destroy the target mechanically.

- destructive effect is governed by **kinetic energy** \(E_k = \tfrac{1}{2} m v^2\),
- penetration is more directly tied to **momentum** \(p = m v\).

Relatively "slow" but heavy projectiles (velocities \(10–10^3\,\text{km/s}\)) often penetrate deep into the target, while ultrafast particle streams are stopped in a thin surface layer.

### 4.1 Limitations of kinetics in space

The main drawback of kinetic projectiles is their **low speed relative to c**, with several consequences:

- long flight time at large distances ⇒ the target has time to maneuver,
- limited effective range (practically \(10^3–10^4\,\text{km}\) for ship‑to‑ship combat),
- low rate of fire of heavy accelerators, limited number of main barrels.

The original text notes that **acceleration tracks** for projectiles are extremely expensive and their length and cost grow roughly with the 2nd–3rd power of the demanded exit velocity (depending on engineering limits).

### 4.2 Basic model of projectile acceleration

For a simple model of a kinetic gun with **constant driving force** \(F\) and projectile mass \(m\), accelerated over distance \(x\), we can write:

\[
 v(x) = \sqrt{\frac{2 F x}{m}}.
\]

The original LaTeX text contains a simplified form (without the factor 2), which we include here only as an illustration:

\[
 v(x) \propto \sqrt{\frac{F x}{m}}.
\]

If the system is limited by **input power** \(P\) instead of force (for example by power limits of the current source in a railgun), different scaling appears – see the notes in `Vykony_kosmickych_lodi.md` and the chapter `07_Physical_model_and_equations.md`.

### 4.3 Projectiles with nuclear warheads

Nuclear warheads are fragile under acceleration – conventional designs do not survive accelerations of order \(10^5 g\), common for kinetic projectiles. This limits the possibility of firing them from classic guns at very short ranges.

An alternative is **compaction upon impact** – a long rod‑shaped projectile made of fissile material with low critical mass and a possible fusion booster. On impact into a massive body (asteroid, planetary crust) the rod can be compressed enough to reach nuclear supercriticality and explode inside the target.

Against lightly built spacecraft this concept faces a problem:

- defensive shields often consist of **thin layers** that tend to vaporize and disperse the projectile rather than compress it coherently.

Nuclear projectiles are therefore more suitable for bombardment of **massive bodies** (military bases in asteroids, etc.), where the depth of the explosion dramatically affects the efficiency of coupling to shock waves.

---

## 5. Guided missiles

Guided missiles are essentially **small spacecraft** with:

- their own propulsion system,
- maneuvering capability,
- usually a nuclear or similarly sophisticated warhead.

### 5.1 The propulsion problem for missiles

For use in space combat it is difficult to build **small but very powerful engines** with high specific impulse:

- most efficient space engines perform better at **larger scales** (especially nuclear engines),
- a missile therefore typically either has high thrust but very low specific impulse (chemical engines),
- or uses exotic, expensive fuel (nuclear isomers) with limited availability.

Examples:

- **Chemical engines** – high thrust, low specific impulse, \(\Delta v \lesssim 40\,\text{km/s}\) at the cost of a very unfavourable fuel‑to‑payload ratio.
- **Engines based on nuclear isomers** – compact nuclear reactors with high specific impulse but extremely expensive fuel preparation.

### 5.2 Cost and role of guided missiles

Guided missiles differ from other weapons in several ways:

- they are **very expensive** – cost comparable to a small spacecraft,
- their effectiveness relies on **sophisticated nuclear warheads**, often with various ways of directing the released energy (EMP, nuclear‑pumped lasers, plasma jets),
- therefore they are used **sparingly** and in carefully prepared tactical situations.

Their advantages include:

- at relatively low accelerations they can carry **complex warheads and electronics**,
- they can execute complex maneuvers and coordinated "wolf‑pack" attacks against ships with limited maneuverability.

---

## 6. Link to tactics and defence

For a more detailed discussion of how each weapon type is used, see:

- [01_Tactics](01_Tactics.md) – long‑range vs. close‑range combat, broadside salvos, preparation for the culminating clash,
- [03_Defense](03_Defense.md) – methods of protection against lasers, particle beams, kinetic projectiles and missiles.

Physical and mathematical models (range, dispersion, hit probability, shield erosion) will be systematically summarized in the chapter [07_Physical_model_and_equations](07_Physical_model_and_equations.md), including tables of parameters and derivations where useful.
