

class MotorType{
	double efficiency;  // how much power is used for acceleration of propelant
	double veMin;       // minimal exhaust velocity 
	double veMax;       // maximal exhaust velocity
	bool exhaustFuel;   // if true the burned fuel is added to propellant mass
	FuelType      fuel;
	PropelantType propelant;	
}

class Motor : public ShipComponent {
	public:
	MotorType typ;
	double size;
	double thrust;
	double power;
	double consumption;
	double mass;

	int    npoint;     // number of anchor points
	int    * points;   // anchor points - refers to points of the Ship
	double * weights;  // distribution of thrust among anchor points -
	Vec3d  * dirs;     // normalized thrust vectors

}

