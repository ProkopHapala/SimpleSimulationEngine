

class GunType{
	double recoil;
	

}

class Gun : public ShipComponent{
	public:
	GunType typ;
	
	int    npoint;     // number of anchor points
	int    * points;   // anchor points - refers to points of the Ship

}
