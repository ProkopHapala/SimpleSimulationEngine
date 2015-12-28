

class ShipComponent{
	double mass;
}

class Tank : public ShipComponent {
	public:
	Comodity * typ;
	double mass;
	double volume;
	double filled;
}

class Pipe : public ShipComponent {
	public:
	ShipComponent * a;
	ShipComponent * b;
}

class Hub : public ShipComponent {
	public:
	int npipes;      // number of pipes going to hub;
	Pipe * pipes;    // 
}



