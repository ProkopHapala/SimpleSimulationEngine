import math

# --- TARGET DEFINITIONS ---
# We define the physical properties for different armor types on different body parts.
# yield_strength_Pa: Pressure to deform/break. Flesh is weak, steel is strong.
# thickness_m: Armor thickness.
# blunt_mass_kg: The effective mass of the body part for calculating blunt trauma.
# damage_multiplier: How lethal a hit to this part is after penetration.
TARGETS = {
    'head': {
        'unarmored': {'yield_strength_Pa': 50_000,      'thickness_m': 0.01,  'blunt_mass_kg': 5.0, 'damage_multiplier': 3.0},
        'cloth':     {'yield_strength_Pa': 500_000,     'thickness_m': 0.01,  'blunt_mass_kg': 5.5, 'damage_multiplier': 3.0},
        'steel':     {'yield_strength_Pa': 500_000_000, 'thickness_m': 0.002, 'blunt_mass_kg': 8.0, 'damage_multiplier': 3.0},
    },
    'torso': {
        'unarmored': {'yield_strength_Pa': 50_000,      'thickness_m': 0.02,   'blunt_mass_kg': 12.0, 'damage_multiplier': 1.5},
        'cloth':     {'yield_strength_Pa': 500_000,     'thickness_m': 0.015,  'blunt_mass_kg': 13.0, 'damage_multiplier': 1.5},
        'steel':     {'yield_strength_Pa': 500_000_000, 'thickness_m': 0.0025, 'blunt_mass_kg': 18.0, 'damage_multiplier': 1.5},
    },
    'limb': {
        'unarmored': {'yield_strength_Pa': 50_000,      'thickness_m': 0.015,  'blunt_mass_kg': 4.0, 'damage_multiplier': 0.8},
        'cloth':     {'yield_strength_Pa': 500_000,     'thickness_m': 0.01,   'blunt_mass_kg': 4.5, 'damage_multiplier': 0.8},
        'steel':     {'yield_strength_Pa': 400_000_000, 'thickness_m': 0.0018, 'blunt_mass_kg': 6.0, 'damage_multiplier': 0.8},
    }
}

def calculate_damage(
    impact_mass_kg: float,
    impact_velocity_mps: float,
    impact_type: str,  # 'cut' or 'pierce'
    target_body_part: str,
    target_armor_type: str,
    edge_length_m: float = 0.0,
    point_angle_deg: float = 0.0
):
    """
    Calculates damage from a single impact based on its physical properties.
    """
    # 1. Calculate Initial Impact Properties
    kinetic_energy = 0.5 * impact_mass_kg * (impact_velocity_mps ** 2)
    
    # 2. Get Target Properties
    try:
        target_stats = TARGETS[target_body_part][target_armor_type]
    except KeyError:
        return {"error": "Invalid target body part or armor type specified."}

    # 3. Calculate Crushing (Blunt) Damage
    # This represents the shockwave/trauma transferred through the armor.
    # It is always applied, regardless of penetration.
    energy_crushing = kinetic_energy * (impact_mass_kg / (impact_mass_kg + target_stats['blunt_mass_kg']))

    # 4. Calculate Armor Penetration and Penetrating Damage
    energy_penetrating = 0.0
    penetration_depth_mm = 0.0

    if impact_type == 'cut':
        if edge_length_m <= 0: return {"error": "Edge length must be positive for a 'cut' impact."}
        
        # Toughness against a cutting attack
        armor_toughness_c = target_stats['yield_strength_Pa'] * (target_stats['thickness_m'] ** 2)
        energy_density_c = kinetic_energy / edge_length_m

        if energy_density_c > armor_toughness_c:
            # Armor is defeated. Calculate energy remaining after the cut.
            excess_density = energy_density_c - armor_toughness_c
            energy_penetrating = excess_density * edge_length_m
            # Simple model for wound depth
            penetration_depth_mm = (energy_penetrating / kinetic_energy) * 50.0 # 50mm is a deep gash

    elif impact_type == 'pierce':
        if point_angle_deg <= 0: return {"error": "Point angle must be positive for a 'pierce' impact."}

        # Model piercing as a series of steps, as the area grows with depth
        remaining_energy = kinetic_energy
        step_m = 0.0001 # 0.1 mm steps
        current_depth_m = 0.0
        
        while remaining_energy > 0 and current_depth_m < target_stats['thickness_m']:
            current_depth_m += step_m
            # Radius of the cone at the current depth
            radius = current_depth_m * math.tan(math.radians(point_angle_deg / 2))
            # Area of the cone's cross-section at this depth
            area = math.pi * (radius ** 2)
            # Energy required to push through this step
            energy_for_step = target_stats['yield_strength_Pa'] * area
            
            if remaining_energy >= energy_for_step:
                remaining_energy -= energy_for_step
            else:
                # Not enough energy to complete the step, stop here
                current_depth_m -= step_m # revert last step
                remaining_energy = 0
        
        penetration_depth_mm = current_depth_m * 1000
        # Penetrating energy is what's left after passing completely through the armor
        if penetration_depth_mm >= target_stats['thickness_m'] * 1000:
            energy_penetrating = remaining_energy

    # 5. Calculate Final Damage Score
    # The total raw damage is the trauma from the impact plus any damage from penetration.
    total_raw_damage = energy_crushing + energy_penetrating
    final_damage_score = total_raw_damage * target_stats['damage_multiplier']
    
    return {
        "Target": f"{target_armor_type.capitalize()} {target_body_part.capitalize()}",
        "Initial KE (J)": round(kinetic_energy, 2),
        "Crushing Damage (J)": round(energy_crushing, 2),
        "Penetrating Damage (J)": round(energy_penetrating, 2),
        "Armor Penetration (mm)": round(penetration_depth_mm, 2),
        "FINAL DAMAGE SCORE": round(final_damage_score, 2)
    }

# --- EXAMPLE SCENARIOS ---

# Define some generic impacts to test the model
impacts = {
    "Light Sword Slash":  {'mass': 0.6, 'vel': 25, 'type': 'cut',    'edge': 0.6,  'angle': 0},
    "Heavy Axe Chop":     {'mass': 1.8, 'vel': 20, 'type': 'cut',    'edge': 0.15, 'angle': 0},
    "Fast Dagger Stab":   {'mass': 0.3, 'vel': 15, 'type': 'pierce', 'edge': 0,    'angle': 15},
    "Spear Power Thrust": {'mass': 2.3, 'vel': 13, 'type': 'pierce', 'edge': 0,    'angle': 10},
}

test_scenarios = [
    ('torso', 'unarmored'),
    ('torso', 'cloth'),
    ('torso', 'steel'),
    ('head', 'steel'),
]

# Run the simulations
for name, impact in impacts.items():
    print(f"\n--- IMPACT: {name.upper()} ---", impact)
    for part, armor in test_scenarios:
        result = calculate_damage(
            impact_mass_kg=impact['mass'],
            impact_velocity_mps=impact['vel'],
            impact_type=impact['type'],
            target_body_part=part,
            target_armor_type=armor,
            edge_length_m=impact['edge'],
            point_angle_deg=impact['angle']
        )
        print(result)