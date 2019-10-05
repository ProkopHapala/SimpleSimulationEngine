

//   From:
//   ftp://ftp.lowell.edu/pub/elgb/astorb.html


#ifndef  Asteroid_h
#define  Asteroid_h

#include <math.h>

constexpr const double DEG_2_RAD = M_PI/180;

//double julian_date_conv( double day, double month, double year ){
double jd_conv( double day, double month, double year ){
    //double m,y,a,b,c,d,jd;
    double m = month;
    double y = year;
    if ( m < 3 ){
        y=y-1;
        m=m+12;
    }
    double a  = floor    ( y * 0.01 );
    double b  = 2-a+floor( a * 0.25 );
    double c  = floor(365.25*y);
    double d  = floor(30.6001*(m+1));
    double jd = b+c+d+day+1720994.5;
    return jd;
}

double true_anomaly( double ma, double ec ){
    // https://en.wikipedia.org/wiki/True_anomaly
    // fast approximation https://en.wikipedia.org/wiki/True_anomaly#From_the_mean_anomaly
    // https://en.wikipedia.org/wiki/Equation_of_the_center#Series_expansion
    double errConv = 1e-6
    double e,diff,f;
    e=ma;
    diff=1.0;
    while(diff>errConv) {
        e = ma + ec*sin(e);
        diff= abs(e-ec*sin(e)-ma);
    }
    f=2*atan(sqrt((1+ec)/(1-ec))*tan(e*0.5));
    return f;
}

void get_position_2d( double a, double e, double i, double lo, double so, double ma ){

    //double inc = i*DEG_2_RAD;
    //double l   = lo*DEG_2_RAD;
    //double s   = so*DEG_2_RAD;
    //double m   = ma*DEG_2_RAD;

    double p=a*(1-e*e);
    double ang=true_anomaly(m,e); // Mean-anomaly to True anomaly
    double cos_ang = cos(ang);
    double sin_ang = sin(ang);
    double rad=p/(1+e*cos_ang);
    double tx=rad*cos_ang; // longer axis
    double ty=rad*sin_ang; // shorter axis

    double cos_l = cos(l);
    double sin_l = sin(l);

    double cos_s = cos(s);
    double sin_s = sin(s);

    //   calculate the rotation matrix
    xx=  cos_s*cos_l - sin_s*sin_l;
    xy= -sin_s*cos_l - cos_s*sin_l;
    yx=  cos_s*sin_l + sin_s*cos_l*cos_inc;
    yy= -sin_s*sin_l + cos_s*cos_l*cos_inc;

    *x=tx*xx+ty*xy;
    *y=tx*yx+ty*yy;
}

int  draw_orbit( double a, double e, double i, double lo, double so ){

    //double inc = i *DEG_2_RAD;
    //double l   = lo*DEG_2_RAD;
    //double s   = so*DEG_2_RAD;
    double p   = a *(1-e*e);    //  semi-latus rectum   https://en.wikipedia.org/wiki/Ellipse#Polar_form_relative_to_focus

    double cos_l = cos(l);
    double sin_l = sin(l);

    double cos_s = cos(s);
    double sin_s = sin(s);

    double xx=  cos_s*cos_l - sin_s*sin_l;
    double xy= -sin_s*cos_l - cos_s*sin_l;
    double yx=  cos_s*sin_l + sin_s*cos_l*cos_inc;
    double yy= -sin_s*sin_l + cos_s*cos_l*cos_inc;

    double  ang=true_anomaly(m,e); // we can save this ?
    double  rad=p/(1+e*cos(ang));
    double  tx=rad*cos(ang);
    double  ty=rad*sin(ang);

    double astep=3/(a*scale);

    for(double ang=0; ang<2*PI; ang+=astep ){
        double cos_ang = cos(ang);
        double sin_ang = sin(ang);
        double r = p/( 1 + e*cos_ang );
        tx=r*cos_ang;
        ty=r*sin_ang;

        x = x0 + scale*( tx*xx + ty*xy );
        y = y0 + scale*( tx*yx + ty*yy );
        z = z0 + scale*( tx*zx + ty*zy );

        //x = sx/2 + scale*( tx*xx + ty*xy );
        //y = sy/2 - scale*( tx*yx + ty*yy );
        //draw_pixel(bm,sx,sy,x,y,0xff00);
    }
}

struct Asteroid{

    double p;
    double e;   // excentricity
    double m0;  // magnitude ?
    double md;


    double a;
    double lo;
    double so;
    double ma;
    /*
    double xx;
    double xy;
    double yx;
    double yy;
    */
    int type;

    // ---- Functions

    void fromString(char* s, double jd ){
        constexpr const double DEG_2_RAD = M_PI/180;
        //double lo,so;
        char blah[256];
        char date[16];
        sscanf(inp,"%104c %s %lf %lf %lf %lf %lf %lf %[^\n]s", blah, date, &ma, &so, &lo, &i, &e, &a, blah );
        fgetc(inp);
        double q=a*(1-e);
        //if(q<min){ // objects which are too far out won't be included
        //double y,m,d,epoch,period,m_now;
        double d = atof(&date[6]);   date[6]=0;
        double m = atof(&date[4]);   date[4]=0;
        double y = atof(date);
        double epoch  = jd_conv(d,m,y);
        double period = pow(a,1.5) * 365.25;
        //a_list[num_asteroids].e=e;
        i  = i*DEG_2_RAD; // inclination
        so = so*DEG_2_RAD;
        lo = lo*DEG_2_RAD;
        md=2*M_PI/period;
        ma=ma*DEG_2_RAD;
        m0=ma + (jd-epoch)*md;

        p=a*(1-e*e);
        if     (q<1  ){ type=2; }
        else if(q<1.3){ type=1; }
        else          { type=0; }
        //}
    }

    void get_position_2d( double *x,double *y ){
        get_position_2d( a, e, i, lo, so, ma, x,y);
    };

};


// loads the asteroid file and calculate the rotation matrix
// objects which are too far out won't be included

/* draws the terrestrial planet */
void draw_solar_system(char *bm,int sx,int sy,double jd,double scale){
double x,y,m_now;

  // draw Sun
  fill_circle(bm,sx,sy,sx/2,sy/2,8,0xffff00);
  // Earth
  draw_orbit(bm,sx,sy,1.0,0.0167,0.00035,357.8,105.245,scale);
  m_now= 23.05797 +  (jd - 2450840.5 ) * 0.985606 ;
  get_position(1.0,0.0167,0.00035,357.8,105.245,m_now,&x,&y);
  x=sx/2+scale*x;
  y=sy/2-scale*y;
  fill_circle(bm,sx,sy,x,y,5,0x00ffff);

  // Test Earth
  /*
    draw_orbit(bm,sx,sy,1.0,0.0167,0.0,0.0,102.252,scale);
    m_now= 100.15 +  (jd - jd_conv(1.5,1,1960) ) * 0.985606 ;
    get_position(1.0,0.0167,0.0,0.0,102.252,m_now,&x,&y);
    x=sx/2+scale*x;
    y=sy/2-scale*y;
    fill_circle(bm,sx,sy,x,y,5,0x00ffff);
  */

  // Mercury
  draw_orbit(bm,sx,sy,0.3870972,0.2056321,7.00507,48.3333,29.1192,scale);
  m_now = 171.74808 +  (jd - 2450840.5 ) *4.092362;
  get_position(0.3870972,0.2056321,7.00507,48.3333,29.1192,m_now,&x,&y);
  x=sx/2+scale*x;
  y=sy/2-scale*y;
  fill_circle(bm,sx,sy,x,y,4,0xffff);

  // Venus
  draw_orbit(bm,sx,sy,0.7233459,0.0067929,3.39460,76.6870,54.871,scale);
  m_now = 1.72101 +  (jd - 2450840.5 ) * 1.602085;
  get_position(0.7233459,0.0067929,3.39460,76.6870,54.871,m_now,&x,&y);
  x=sx/2+scale*x;
  y=sy/2-scale*y;
  fill_circle(bm,sx,sy,x,y,4,0xffff);

  // Mars
  draw_orbit(bm,sx,sy,1.5237118,0.0934472,1.84990,49.5629,286.4653,scale);
  m_now = 10.2393 +  (jd - 2450840.5 ) * 0.5240224;
  get_position(1.5237118,0.0934472,1.84990,49.5629,286.4653,m_now,&x,&y);
  x=sx/2+scale*x;
  y=sy/2-scale*y;
  fill_circle(bm,sx,sy,x,y,4,0xffff);

}


#endif
