
#ifndef FieldPatch_h
#define FieldPatch_h

#include "Vec3.h"

inline float trinagle_area( float x12, float y12, float x13, float y13 ) {  return x12*y13 - x13*y12; };

class Quad{
	public:
	double value;
	union{
		struct{
			double x1,y1,z1;
			double x2,y2,z2;
			double x3,y3,z3;
			double x4,y4,z4;
		};
		struct{ Vec3d p1,p2,p3,p4; }; // order does matter for default initialized
	};

	// ==== function declarations

	void draw() const;
	Quad(  double value_, const Vec3d& p1_, const Vec3d& p2_, const Vec3d& p3_, const Vec3d& p4_ );

	// ==== inline functions

	inline double area() const {
		double x12 = x1 - x2; double y12 = y1 - y2;
		double x13 = x1 - x3; double y13 = y1 - y3;
		double x42 = x4 - x2; double y42 = y4 - y2;
		double x43 = x4 - x3; double y43 = y4 - y3;
		return fabs( trinagle_area( x12, y12, x13, y13 ) ) + fabs(trinagle_area( x42, y42, x43, y43 ) );
	};
};



class FieldPatch {
	public:
	static constexpr float thresh[4] = { 0.6 };
	int   startkill  = 4;
	float dval       = 0.2;
	float dmax       = 0.3;
	float minarea    = 1000;
	float dheight    = 100;
	float dmheight   = 50;

	int quadCount;

	// ==== function declarations

	void divide_1( int nlevels, int level, const Quad& rc   );
	void divide  ( int nlevels, int level, const Quad& quad );
	int  makeList( int nlevels, const Quad& quad );

};

#endif
