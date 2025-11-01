# `asteroidEngineering.h`

This header file lays the groundwork for mechanics related to asteroid utilization, including their composition, and the physics of moving them. It hints at potential gameplay features involving mining, resource management, and the strategic use of asteroids.

## Data Structures

### `MineralType` Struct

Defines a mineral based on its elemental composition.

*   `std::string name`: The name of the mineral.
*   `std::unordered_map<int, double> elements`: A map where the key is the atomic number of an element and the value is its fractional abundance in the mineral.

### `RockType` Struct

Defines a type of rock, which is a composite of different minerals.

*   `std::string name`: The name of the rock type.
*   `std::unordered_map<MineralType*, double> minerals`: A map defining the mineral composition of the rock.
*   `double hardness`: A value representing how difficult the rock is to cut or drill.
*   `double ablationEnergy`: The energy required to ablate (e.g., with a laser) or evaporate the rock.
*   `SpectralProperties spectal`: The spectral properties of the rock, used for remote sensing.

### `Asteroide` Class

Represents an asteroid as a game object, combining its physical and orbital properties with its resource composition.

*   `Orbit* orbit`: A pointer to the asteroid's `Orbit` object.
*   `double mass`: The total mass of the asteroid.
*   `std::unordered_map<RockType*, Deposit> rocks`: A map detailing the composition of the asteroid in terms of different rock types.

## Functions

### `manuever_planeChange`

`double manuever_planeChange(double angle, double mass, Orbit* orbit, double& mProp, double& deltaV)`

This function calculates the requirements for changing the orbital plane of an object (like an asteroid).

*   **Parameters**:
    *   `double angle`: The desired change in orbital inclination (in radians).
    *   `double mass`: The total mass of the object to be moved.
    *   `Orbit* orbit`: The current orbit of the object.
    *   `double& mProp`: A reference to store the calculated required propellant mass.
    *   `double& deltaV`: A reference to store the calculated required change in velocity (delta-v).
*   **Returns**: The total energy `E` required for the maneuver.
*   **Functionality**: It calculates the necessary delta-v for the plane change, and then uses the Tsiolkovsky rocket equation to determine the propellant mass required, assuming a specific engine performance (exhaust velocity).
