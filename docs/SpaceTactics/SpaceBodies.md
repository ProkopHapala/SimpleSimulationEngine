# `SpaceBodies.h`

This header file defines the fundamental data structures for all physical objects in the simulation, including their orbital properties and physical attributes. It provides the building blocks for representing planets, ships, and asteroids.

## `OrbitalElements` Struct

Defines a celestial body's orbit using the classical set of six orbital elements.

*   `double semi_major`: The semi-major axis of the orbit.
*   `double eccentricity`: The eccentricity of the orbit (0 for circular, <1 for elliptical).
*   `double inclination`: The inclination of the orbital plane.
*   `double arg_periapsis`: The argument of periapsis.
*   `double node_longitude`: The longitude of the ascending node.
*   `double mean_anomaly`: The mean anomaly at a specific epoch.
*   `double period`: The orbital period.
*   `double epoch`: The reference time at which the orbital elements are defined.

## `Orbit` Class

Represents the orbit of a celestial body, providing methods to calculate its position and velocity.

*   `Mat3d rot`: A rotation matrix that transforms from the orbital plane to the world coordinate system.
*   `double L`: The specific angular momentum of the orbit.
*   `double E`: The specific orbital energy.
*   `double eccentricity`, `double semi_major`: Key orbital shape parameters.
*   `double mean_anomaly`, `double period`, `double epoch`: Parameters defining the object's position in the orbit over time.

### Key Methods

*   `Orbit(const OrbitalElements& el)`: Constructor to create an `Orbit` object from a set of `OrbitalElements`.
*   `double speed(double r, double M)`: Calculates the orbital speed at a given radius `r` from the central body of mass `M` using the vis-viva equation.
*   `Vec3d pointAtEpoch(double at_epoch) const`: Calculates the 3D position of the orbiting body at a specific time (`at_epoch`).
*   `void toPoints(...)`: Generates a series of points along the orbit path.

## `SpaceBody` Class

Represents any physical object in the simulation world, such as a planet, ship, or asteroid. It inherits from `RigidBody`.

### Public Members

*   `std::string name`: The name of the body.
*   `double radius`: The physical radius of the body.
*   `Vec3d* trjPos`: A pointer to an array of `Vec3d` storing the pre-calculated positions of the body over time.
*   `Vec3d* trjThrust`: A pointer to an array of `Vec3d` storing the thrust vectors applied to the body at each step of its trajectory.
*   `Orbit* orbit`: A pointer to an `Orbit` object if the body follows a Keplerian orbit.
*   `SpaceBody* orbCenter`: A pointer to the central body around which this body orbits.

### Key Methods

*   `SpaceBody(...)`: Constructors for creating a `SpaceBody`.
*   `Vec3d pointAtEpoch(double epoch) const`: Returns the position of the body at a specific time, calculated from its orbit if available.
*   `Vec3d getTrjPos(int iTrj, double du) const`: Interpolates and returns the position of the body from its pre-calculated trajectory at a fractional index.
*   `Vec3d getThrust(int itrj, double du)`: Interpolates and returns the thrust vector from the `trjThrust` array.
*   `void fromString_astorb(char* s)`: Parses a string from the `astorb` database format to initialize the properties of an asteroid, including its orbital elements.

## `SpaceBodyIntegrator` Class

An ODE solver specifically designed to integrate the trajectory of a single `SpaceBody` under the gravitational influence of a set of central bodies.

*   `SpaceBody* o`: The object whose trajectory is being integrated.
*   `SpaceBody** centers`: An array of pointers to the central bodies providing the gravitational field.
*   `virtual void getDerivODE(...)`: The derivative function that calculates the gravitational forces on the body `o`.
*   `void evalTrj(...)`: Runs the integration to compute the trajectory.
