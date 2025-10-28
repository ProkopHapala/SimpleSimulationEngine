import math

# --- BIOMECHANICAL PARAMETERS (Tweak these to change soldier stats) ---

# Parameters for a two-handed swing
SOLDIER_SWING_PARAMS = {
    "torque_Nm": 250,        # Rotational force in Newton-meters
    "angle_rad": 1.5708      # Swing angle (90 degrees = ~1.5708 radians)
}

# Parameters for a one-handed thrust (fast jab, e.g., with a shield)
SOLDIER_THRUST_1H_PARAMS = {
    "force_N": 150,          # Linear force from arm/shoulder
    "distance_m": 0.8,       # Acceleration distance
    "m_dead_human_kg": 2     # "Wasted" mass of the accelerating arm
}

# Parameters for a two-handed power thrust (full body lunge)
SOLDIER_THRUST_2H_PARAMS = {
    "force_N": 500,          # Linear force from whole body
    "distance_m": 1.5,       # Acceleration distance of the lunge
    "m_dead_human_kg": 15    # "Wasted" mass of the lunging body
}

# --- WEAPON DEFINITIONS ---
# Each weapon is a dictionary with its physical properties.
# L: Total Length (m)
# m_total: Total Mass (kg)
# m_head: Mass of the metal head (kg)
# weapons = [
#     # name, L, m_total, m_head
#     {"name": "Spear",   "L": 2.5, "m_total": 2.0, "m_head": 0.3},
#     {"name": "Pike",    "L": 4.5, "m_total": 5.0, "m_head": 0.2},
#     {"name": "Glaive",  "L": 2.2, "m_total": 2.5, "m_head": 1.0},
#     {"name": "Halberd", "L": 2.0, "m_total": 3.0, "m_head": 1.5},
#     {"name": "Poleaxe", "L": 1.8, "m_total": 3.5, "m_head": 2.2},
# ]

weapons = [
    # name, L, m_head. We'll add shaft_kg_per_m to each.
    {"name": "Spear",   "L": 2.5, "m_head": 0.3, "shaft_kg_per_m": 0.8},
    {"name": "Pike",    "L": 5.0, "m_head": 0.2, "shaft_kg_per_m": 1.2},
    {"name": "Glaive",  "L": 2.0, "m_head": 1.5, "shaft_kg_per_m": 1.0},
    {"name": "Halberd", "L": 2.5, "m_head": 1.2, "shaft_kg_per_m": 1.0},
    {"name": "Poleaxe", "L": 1.8, "m_head": 2.2, "shaft_kg_per_m": 1.2},
]

# --- CALCULATION LOGIC ---

results = []
for weapon in weapons:
    # --- Basic Properties ---
    name = weapon["name"]
    L = weapon["L"]
    # m_total = weapon["m_total"]
    # m_head = weapon["m_head"]
    # m_shaft = m_total - m_head

    m_head = weapon["m_head"]
    m_shaft = L * weapon["shaft_kg_per_m"]
    m_total = m_shaft + m_head

    # --- SWING CALCULATIONS ---
    # Moment of Inertia (I): Resistance to rotation.
    # Formula: I = (1/3 * m_shaft * L^2) + (m_head * L^2)
    I = ((1/3 * m_shaft) + m_head) * (L ** 2)
    
    # Work done by the soldier during the swing
    work_swing = SOLDIER_SWING_PARAMS["torque_Nm"] * SOLDIER_SWING_PARAMS["angle_rad"]
    
    # Final angular velocity (omega) in radians/sec
    # KE_rot = 0.5 * I * omega^2 => omega = sqrt(2 * Work / I)
    omega = math.sqrt(2 * work_swing / I)
    
    # Time to complete the swing
    # alpha = tau / I; theta = 0.5 * alpha * t^2 => t = sqrt(2 * theta / alpha)
    time_swing = math.sqrt(2 * SOLDIER_SWING_PARAMS["angle_rad"] * I / SOLDIER_SWING_PARAMS["torque_Nm"])

    # Velocity of the weapon's head at impact
    v_head = omega * L
    
    # EFFECTIVE ENERGY & MOMENTUM for swing (damage is from the head)
    ke_swing = 0.5 * m_head * (v_head ** 2)
    p_swing = m_head * v_head

    # --- 1H THRUST CALCULATIONS ---
    f_1h = SOLDIER_THRUST_1H_PARAMS["force_N"]
    d_1h = SOLDIER_THRUST_1H_PARAMS["distance_m"]
    m_dead_1h = SOLDIER_THRUST_1H_PARAMS["m_dead_human_kg"]
    
    # Total mass the soldier must accelerate (weapon + his own arm)
    m_accelerated_1h = m_total + m_dead_1h
    work_1h = f_1h * d_1h
    
    # Final velocity of the thrust
    # Work = KE_total = 0.5 * m_accelerated * v^2 => v = sqrt(2 * Work / m_accelerated)
    v_1h = math.sqrt(2 * work_1h / m_accelerated_1h) if m_accelerated_1h > 0 else 0
    
    # Time to complete the thrust
    # d = 0.5 * a * t^2; a = F/m => t = sqrt(2 * d * m / F)
    time_1h = math.sqrt(2 * d_1h * m_accelerated_1h / f_1h) if f_1h > 0 else float('inf')

    # EFFECTIVE ENERGY & MOMENTUM for thrust (damage is from the whole weapon)
    ke_1h = 0.5 * m_total * (v_1h ** 2)
    p_1h = m_total * v_1h

    # --- 2H THRUST CALCULATIONS ---
    f_2h = SOLDIER_THRUST_2H_PARAMS["force_N"]
    d_2h = SOLDIER_THRUST_2H_PARAMS["distance_m"]
    m_dead_2h = SOLDIER_THRUST_2H_PARAMS["m_dead_human_kg"]
    
    m_accelerated_2h = m_total + m_dead_2h
    work_2h = f_2h * d_2h
    v_2h = math.sqrt(2 * work_2h / m_accelerated_2h) if m_accelerated_2h > 0 else 0
    time_2h = math.sqrt(2 * d_2h * m_accelerated_2h / f_2h) if f_2h > 0 else float('inf')
    
    # EFFECTIVE ENERGY & MOMENTUM for 2H thrust
    ke_2h = 0.5 * m_total * (v_2h ** 2)
    p_2h = m_total * v_2h

    time_factor = 2.0;
    # Store all calculated values
    results.append({
        "Name": name,  
        "Total Mass (kg)": m_total,
        "Time Swing (s)":       time_swing*time_factor,
        "KE Swing (J)":         ke_swing,
        "P Swing (kg·m/s)":     p_swing,
        "Time 1H Thrust (s)":   time_1h*time_factor,
        "KE 1H Thrust (J)":     ke_1h,
        "P 1H Thrust (kg·m/s)": p_1h,
        "Time 2H Thrust (s)":   time_2h*time_factor,
        "KE 2H Thrust (J)":     ke_2h,
        "P 2H Thrust (kg·m/s)": p_2h,
    })

# --- DISPLAY RESULTS ---

# def print_results_table(results_data):
#     """Prints a formatted table of the results."""
#     # Define headers and their widths
#     headers = {
#         "Name": 10, 
#         "T_Sw": 6, "KE_Sw": 7, "P_Sw": 7, "|": 2,
#         "T_1H": 6, "KE_1H": 7, "P_1H": 7, "|": 2,
#         "T_2H": 6, "KE_2H": 7, "P_2H": 7
#     }
#     header_str = "".join([f"{h:<{w}}" for h, w in headers.items()])
#     print(header_str)
#     print("-" * len(header_str))

#     # Map keys from results dictionary to the short headers
#     key_map = {
#         "Name": "Name", 
#         "T_Sw": "Time Swing (s)",     "KE_Sw": "KE Swing (J)",     "P_Sw": "P Swing (kg·m/s)",
#         "T_1H": "Time 1H Thrust (s)", "KE_1H": "KE 1H Thrust (J)", "P_1H": "P 1H Thrust (kg·m/s)",
#         "T_2H": "Time 2H Thrust (s)", "KE_2H": "KE 2H Thrust (J)", "P_2H": "P 2H Thrust (kg·m/s)"
#     }

#     for res in results_data:
#         row_str = ""
#         for h, w in headers.items():
#             if h == "|":
#                 row_str += f"{'|':<{w}}"
#                 continue
#             key = key_map[h]
#             value = res[key]
#             if isinstance(value, float):
#                 row_str += f"{value:<{w}.2f}"
#             else:
#                 row_str += f"{value:<{w}}"
#         print(row_str)

# print("Polearm Performance Model (T=Time, KE=Kinetic Energy, P=Momentum)")
# print("       |  Swing               | 1H Thrust         | 2H Thrust")
# print_results_table(results)


# --- DISPLAY RESULTS ---
def print_results_table(results_data):
    headers = {
        "Name": 10, "M Tot": 7, 
        "T_Sw": 6, "KE_Sw": 7, "P_Sw": 7, "|": 2,
        "T_1H": 6, "KE_1H": 7, "P_1H": 7, "|": 2,
        "T_2H": 6, "KE_2H": 7, "P_2H": 7
    }
    header_str = "".join([f"{h:<{w}}" for h, w in headers.items()])
    print(header_str)
    print("-" * len(header_str))
    key_map = {
        "Name": "Name", "M Tot": "Total Mass (kg)",
        "T_Sw": "Time Swing (s)",     "KE_Sw": "KE Swing (J)",     "P_Sw": "P Swing (kg·m/s)",
        "T_1H": "Time 1H Thrust (s)", "KE_1H": "KE 1H Thrust (J)", "P_1H": "P 1H Thrust (kg·m/s)",
        "T_2H": "Time 2H Thrust (s)", "KE_2H": "KE 2H Thrust (J)", "P_2H": "P 2H Thrust (kg·m/s)"
    }
    for res in results_data:
        row_str = ""
        for h, w in headers.items():
            if h == "|": row_str += f"{'|':<{w}}"; continue
            key = key_map[h]; value = res[key]
            row_str += f"{value:<{w}.2f}" if isinstance(value, float) else f"{value:<{w}}"
        print(row_str)

print("Polearm Performance Model (Modular Definition) (T=Time, KE=Kinetic Energy, P=Momentum)")
print("                | Swing              | 1H Thrust (Jab)    | 2H Thrust (Power)")
print_results_table(results)