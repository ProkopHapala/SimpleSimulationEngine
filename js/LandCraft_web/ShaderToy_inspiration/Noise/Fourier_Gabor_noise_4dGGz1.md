

## Buffer A

```glsl
// do your operation in spectral domain here. 
// Or tune #PROF #PHASE line 35-36

#define SIZE (iResolution.x/2.-30.) //Size must be changed in each tab.

bool keyPress(int ascii) { return (texture(iChannel2,vec2((.5+float(ascii))/256.,0.25)).x > 0.); }
float rand(vec2 uv) { return fract(1e5*sin(dot(uv,vec2(17.4,123.7)))); }     // point -> rnd [0,1[
#define srnd(x) (2.*rand(x)-1.)
float gauss(float x) { return exp(-.5*x*x); }
#define ang(a)  vec2(cos(a), sin(a))                         // used to set the complex exp(i.a)
vec2 cmul (vec2 a,float b) { return mat2(a,-a.y,a.x) * vec2(cos(b),sin(b)); } // complex a*exp(i.b)


void mainImage( out vec4 O, vec2 U )
{
    O*=0.;
    vec2 R = iResolution.xy;
    if ( U==vec2(.5)) {
        if (iFrame==0) O.zw = vec2(0,0);
        else           O.zw = texture(iChannel1, U/R).zw;  
        if ( keyPress(32) ) 
            if (iMouse.x/R.x<.5) O.z = mod(O.z+.1, 5.) ; // persistant key flag for left window
            else                 O.w = mod(O.w+.1, 5.) ; // persistant key flag for right window
        return;
    }
    
    U -= .5;  // freq 0 must exist
    
    U = 2.*U-SIZE;
    vec2 X = U/SIZE, T,
         M = 2.*iMouse.xy/R-1.;
    float I=1., l = length(X), F, s = sign(-X.x); // s to help making a symmetric spectrum (phases included !)
    
    
    // --- your custom Fourier-space function here ------------
#define PROF  32 // (32) spectrum profile ( = Fourier modulus)
#define PHASE 30 // (30) spectrum phases
    
#if 0 // 0:  scale shape at zoom    1: scale fourier at zoom
    X *= SIZE/256.;  l *= SIZE/256., I *= SIZE/256.; 
#endif    
    
                                     // --- modulus profile here
#if PROF==0
    F = 1.;                                     // flat
#elif PROF==1
    F = gauss(l/.05)*10.;                       // gauss
#elif PROF==101
    F = gauss(l/.01*M.x)*10.;                   // tunable gauss
#elif PROF==11
    F = exp(-l/.05)*10.;                        // exp
#elif PROF==111
    F = exp(-l/.02*(1.+M.x))*10.;               // tunable exp
#elif PROF==2
    float l1 = length(X-vec2(.07,0)),
          l2 = length(X+vec2(.07,0));  
    F = ( gauss(l1/.02)+gauss(l2/.02) )*10.;    // bi-lobe
  //l1 = length(X-vec2(.1,.05)),
  //l2 = length(X+vec2(.1,.05)); 
  //F += ( gauss(l1/.015)+gauss(l2/.015) )*5.;  // additionnal bi-lobe
#elif PROF==30
    F = step(.5,l); // /(1e-5+l*l);             // --- white-disc (blue noise)
#elif PROF==3
    F = gauss(abs(l-.12)/.005)*10.;             // --- ring (blue noise)
#elif PROF==31
 // F = gauss(abs(l-.22)/.002)*10.;             // variant
 // F = gauss(abs(l-.5)/.05)*2.;                // variant: 1/2 thinnest blue noise
    F = gauss(abs(l-.95)/.01)*4.;               // variant: thinnest blue noise
#elif PROF==32
    F = gauss(abs(l-.32-.1*cos(.1*iTime))/.002)*10.;             // variant
#elif PROF==4
    F = gauss(abs(l-.12)/.007)*10.*gauss(length(X*vec2(.1,1))/.03)*3.;  
#elif PROF==5
    F = fract(sin(dot(U, vec2(12.9898, 78.233)))* 43758.5453) *2.-1.; // random
  //F = fract(sin(dot(U, vec2(13, 79)))* 4e5) *2.-1.; 
#elif PROF==6                                   // --- let generate aliasing ! :-p
  //F = gauss(abs(l-.95)/.002)*10.; // for qualibration.
    l = length(X-.8); l *= SIZE/256.;
    F += gauss(l/.05)*10.;             
    l = length(X+.8); l *= SIZE/256.;
    F += gauss(l/.05)*10.;             
#elif PROF==61                                  // variant
    F = smoothstep(1.,1.02,l-.0)*100.;    // try -.1, -.35 
#elif PROF==70                                  // structure phases
    F = smoothstep(.01,0.,abs(abs(X.x)-.3)) * gauss(X.y/.1) * 3.;
#endif
 
                                     // --- phases here ( 0 for direct Fourier transform)
  //vec2 P = ang(6.2832*rand(U));               // default: random phases
    vec2 P = ang(6.2832*rand(X*s)*s);           // with phase symmetry
#if PHASE==0
    T = vec2(1,0);                              // no phases ( all 0 )
#elif PHASE==1
    T = P;                                      // random phases
#elif PHASE==20
    T = ang(6.2832*length(30.*X));              // correlated phase: linear 
#elif PHASE==21
    T = ang(6.2832*length(sin(30.*X)));         // correlated phase
#elif PHASE==212
    T = ang(6.2832*length(sin(100.*M*X)));      // variant with mouse gain
#elif PHASE==22
    T = normalize( vec2(abs(X.x),X.y)-vec2(.07,0))*.5; // rotating around (.07,0) - use with PROF=2
#elif PHASE==3
    T = cmul(P,2.*iTime*s);                     // phase shift with time (X biased)   
#elif PHASE==30
    float t = .32*iTime, a = fract(t);          // phase shift with time (morph)
    T = ang(6.2832*mix(srnd(X*s+floor(t)),
                       srnd(X*s+ceil (t)),
                       a)  *s / sqrt(a*a+(1.-a)*(1.-a))); // preserve variance along time
#elif PHASE==31
    T = cmul(P,2.*iTime*s*sqrt(30.*abs(X.x))); // dispersive phase shift 
  //T = cmul(P,2.*iTime*s*sqrt(1./(1e-5+abs(X.x)))); 
#endif
    
    
    O = vec4(T*F,0,0)*sqrt(I); //  *SIZE;  
}
```

## Buffer B

```glsl
// invFourier transform 

// Horizontal + Vertical Discrete Fourier Transform of the input 
// 2 passes pipelined : in -> buf.zw -> buf.xy -> out
// ( adapted from  Flyguy's https://www.shadertoy.com/view/MscGWS# )


#define SIZE (iResolution.x/2.-30.) //Size must be changed in each tab.

//#define tex(ch,x,y) texture(ch, vec2(x,y)/iResolution.xy )
#define tex(ch,x,y)  texelFetch(ch, ivec2(x,y), 0)

vec2 cmul (vec2 a,float b) { return mat2(a,-a.y,a.x)*vec2(cos(b), sin(b)); } // complex a*exp(i.b)
#define W(uv)   mod(uv+SIZE/2.,SIZE)                    // wrap [-1/2,1/2] to [0,1]


void mainImage( out vec4 O, vec2 U )
{
    O-=O; 
    
    if(U.x > SIZE || U.y > SIZE) return;

    for(float n = 0.; n < 1000.; n++)  {
        if (n>=SIZE) break;
        float m = W(n);       // W to warp 0,0 to mid-window.
        vec2 xn = tex(iChannel0, m+.5, U.y).xy,
             yn = tex(iChannel1, U.x, m+.5).zw,
             a =  6.2831853 *  W(U-.5) * n/SIZE;
        
        O.zw += cmul(xn, a.x);
        O.xy += cmul(yn, a.y);
    }
    
    O.zw /= SIZE;
}
```

## Image

```glsl
// application of https://www.shadertoy.com/view/4s3GDs

// set you module and phase in Buf A

#define SIZE (iResolution.x/2.-30.) //Size must be changed in each tab.

//Display modes.     Tuned by pressing SPACE after clicking left or right window
#define MAGNITUDE 0.
#define PHASE     1.
#define COMPONENT 2.
#define REAL      3.
#define IMAG      4.

//Scaling
#define LOG 0
#define LINEAR 1

#define MAG_SCALE LINEAR

bool keyToggle(int ascii) { return (texture(iChannel1,vec2((.5+float(ascii))/256.,0.75)).x > 0.);}

vec4 rainbow(float x)  { return .5 + .5 * cos(6.2832*(x - vec4(0,1,2,0)/3.)); }
vec4 rainbow(vec2 C)   { return rainbow(atan(C.y,C.x)/3.1416 + .5); }

vec4 paintDFT(vec2 F, float mode) {
    // F /= SIZE;
    return 
         mode == MAGNITUDE 
     #if   MAG_SCALE == LOG
                           ?  vec4(log(length(F)))
     #elif MAG_SCALE == LINEAR
                           ?  vec4(length(F))
     #endif
       : mode == PHASE     ?  rainbow(F)        
       : mode == COMPONENT ?  .5+.5*vec4(F, 0,0)
       : mode == REAL      ?  .5+.5*vec4(F.x)
       : mode == IMAG      ?  .5+.5*vec4(F.y)
       : vec4(-1); // error
}


void mainImage( out vec4 O,  vec2 uv )
{    
    vec2 R = iResolution.xy;
    
    vec2 pixel = ( uv - R/2.) / SIZE  + vec2(2,1)/2.,
         tile  = floor(pixel),
         stile = floor(mod(2.*pixel,2.));    
	     uv = fract(pixel) * SIZE / R;
    O-=O;

    vec2 DISPLAY_MODE = floor(texture(iChannel3, .5/R).zw); // persistant key flag.
    if (tile.y==-1. && abs(tile.x-.5)<1.) {   // buttons displaying current flags value
        for (float i=0.; i<5.; i++) 
            O += smoothstep(.005,.0,abs(length(uv*R/SIZE-vec2(.2+i/7.,.97))-.025));
        float v = tile.x==0. ? DISPLAY_MODE[0] : DISPLAY_MODE[1];
        O.b += smoothstep(.03,.02,length(uv*R/SIZE-vec2(.2+v/7.,.97)));
    }
    
    if (keyToggle(64+6)) // 'F' 
        O += paintDFT(texture(iChannel2, fract(uv)).xy, DISPLAY_MODE[1]); // tiled display
    else {  
        if(tile == vec2(0,0))  // Input spectrum (Left)
            O += paintDFT(texture(iChannel3, uv).xy, DISPLAY_MODE[0]);

        if(tile == vec2(1,0))  // Output DFT (Right)
            O += paintDFT(texture(iChannel2, uv).xy, DISPLAY_MODE[1]);
        }
}
```
