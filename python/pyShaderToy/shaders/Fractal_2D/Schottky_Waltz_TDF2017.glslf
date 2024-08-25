// Copy From Here
/*
https://www.shadertoy.com/view/XslyzH
[TDF2017] Schottky Waltz 

 Tags: 2d, fractal, kleiniangroups, schottky
Created by soma_arc in 2017-02-19

Tokyo Demo Fest 2017 GLSL Graphics compo 1st place.
Play with circles.

*/


/*
Created by soma_arc - 2017
This work is licensed under Creative Commons Attribution-ShareAlike 3.0 Unported.
*/

// from Syntopia http://blog.hvidtfeldts.net/index.php/2015/01/path-tracing-3d-fractals/
vec2 rand2n(vec2 co, float sampleIndex) {
    vec2 seed = co * (sampleIndex + 1.0);
	seed+=vec2(-1,1);
    // implementation based on: lumina.sourceforge.net/Tutorials/Noise.html
    return vec2(fract(sin(dot(seed.xy ,vec2(12.9898,78.233))) * 43758.5453),
                fract(cos(dot(seed.xy ,vec2(4.898,7.23))) * 23421.631));
}

vec2 circleInvert(vec2 pos, vec3 circle){
	return ((pos - circle.xy) * circle.z * circle.z)/(length(pos - circle.xy) * length(pos - circle.xy) ) + circle.xy;
}

mat2 getRotationMat2(float angleRadians){
	return mat2(cos(angleRadians), -sin(angleRadians),
                sin(angleRadians), cos(angleRadians));
}
const mat2 UNI_MAT = mat2(1, 0, 0, 1);

const float PI = 3.1415926535;
const float PI_2 = 3.1415926535/2.;
const float PI_4 = 3.1415926535/4.;


const float R = 100. * sqrt(2.) / 2.;
vec3 c1 = vec3(100, 0, 0);
vec3 c3 = vec3(-100, 0, 0);
vec3 c2 = vec3(0, 100, 0);
vec3 c4 = vec3(0, -100, 0);
vec3 c5 = vec3(0);

float l1Angle = PI_4;
vec2 l1p = vec2(500, 0);
mat2 l1m;
mat2 l1mInv;

float l2Angle = PI + PI_4;
vec2 l2p = vec2(-800, 0);
mat2 l2m;
mat2 l2mInv;

float l3Angle = PI + PI_4;
vec2 l3p = vec2(-1500, 0);
mat2 l3m;
mat2 l3mInv;

float l4Angle = PI_2 + PI_4;
vec2 l4p = vec2(0, -1500);
mat2 l4m;
mat2 l4mInv;

bool enableL1 = false;
bool enableL2 = false;
bool enableL3 = false;
bool enableL4 = false;

const int MAX_ITERATIONS = 35;
float IIS(vec2 pos){
    float loopNum = 0.;
	bool cont = false;
	for(int i = 0 ; i < MAX_ITERATIONS ; i++){
		cont = false;
		
        if(length(pos - c1.xy) < c1.z){
			pos = circleInvert(pos, c1);
			cont = true;
            loopNum++;
		}else if(length(pos - c2.xy) < c2.z){
			pos = circleInvert(pos, c2);
			cont = true;
            loopNum++;
		}else if(length(pos - c3.xy) < c3.z){
			pos = circleInvert(pos, c3);
			cont = true;
            loopNum++;
		}else if(length(pos - c4.xy) < c4.z){
			pos = circleInvert(pos, c4);
			cont = true;
            loopNum++;
		}else if(length(pos - c5.xy) < c5.z){
			pos = circleInvert(pos, c5);
			cont = true;
            loopNum++;
		}
        
        if(enableL1){
        	pos -= l1p;
        	pos = l1mInv * pos;
        	if(pos.x > 0.){
        		pos.x *= -1.;
            	loopNum++;
            	cont = true;
        	}
        	pos = l1m * pos;
        	pos += l1p;
        }
        
        if(enableL2){
        pos -= l2p;
        pos = l2mInv * pos;
        if(pos.x > 0.){
        	pos.x *= -1.;
            loopNum++;
            cont = true;
        }
        pos = l2m * pos;
        pos += l2p;
        }
        
        if(enableL3){
        pos -= l3p;
        pos = l3mInv * pos;
        if(pos.x > 0.){
        	pos.x *= -1.;
            loopNum++;
            cont = true;
        }
        pos = l3m * pos;
        pos += l3p;
        }
         
        if(enableL4){
        	pos -= l4p;
        	pos = l4mInv * pos;
        	if(pos.x > 0.){
        		pos.x *= -1.;
            	loopNum++;
				cont = true;
        	}
        	pos = l4m * pos;
        	pos += l4p;
        }
        
		if(cont == false) break;
	}

	return loopNum;
}

vec3 hsv2rgb(vec3 c)
{
    vec4 K = vec4(1.0, 2.0 / 3.0, 1.0 / 3.0, 3.0);
    vec3 p = abs(fract(c.xxx + K.xyz) * 6.0 - K.www);
    return c.z * mix(K.xxx, clamp(p - K.xxx, 0.0, 1.0), c.y);
}

//w: start time
//s: duration
float scene(in float t, in float w, in float s){
    return clamp(t - w, 0.0, s) / s;  
}


float expEasingIn(float t){
    return pow( 2., 13. * (t - 1.) );
}
float expEasingOut(float t) {
	return -pow( 2., -10. * t) + 1.;
}

float circEasingInOut(float t){
	t /= .5;
	if (t < 1.) return -.5 * (sqrt(1. - t*t) - 1.);
	t -= 2.;
	return .5 * (sqrt(1. - t*t) + 1.);
}

float circEasingIn(float t){
	return -  (sqrt(1. - t*t) - 1.);
}

const float DISPLAY_GAMMA_COEFF = 1. / 2.2;
vec3 gammaCorrect(vec3 rgb) {
  return vec3((min(pow(rgb.r, DISPLAY_GAMMA_COEFF), 1.)),
              (min(pow(rgb.g, DISPLAY_GAMMA_COEFF), 1.)),
              (min(pow(rgb.b, DISPLAY_GAMMA_COEFF), 1.)));
}

const float SAMPLE_NUM = 20.;
void mainImage(out vec4 fragColor, in vec2 fragCoord){
    vec3 sum = vec3(0);
    float t = mod(iTime, 20.);
	float ratio = iResolution.x / iResolution.y / 2.0;
    
    enableL1 = enableL2 = enableL3 = enableL4 = false;

	float start = 0.;
    c1.z = mix(0., R, expEasingIn(scene(t, start, 1.)));
    c2.z = mix(0., R, expEasingIn(scene(t, start + .1, 1.)));
	c3.z = mix(0., R, expEasingIn(scene(t, start + .3, 1.)));
	c4.z = mix(0., R, expEasingIn(scene(t, start + .5, 1.)));
    
    float rotationStart = start + 1.7;
    float scaleFactor = 150.;    
    float theta = mix(0., 2. * PI + PI_4, circEasingInOut(scene(t, rotationStart, 1.)));
    scaleFactor += mix(0., 800., circEasingInOut(scene(t, rotationStart, 1.)));
        
	l1p.x = mix(1500., 200., expEasingIn(scene(t, rotationStart, 1.)));
	l2p.x = mix(-1500., -200., expEasingIn(scene(t, rotationStart, 1.5)));
    enableL1 = t > rotationStart + .5;
    enableL2 = t > rotationStart + .5;
    
	float lineStart = rotationStart + 1.5;
    l1Angle += mix(0., PI_2 , circEasingInOut(scene(t, lineStart, 1.)));
    l2Angle += mix(0., PI_2 , circEasingInOut(scene(t, lineStart, 1.)));
    float lpRotateAngle1 = mix(0., PI_2, circEasingInOut(scene(t, lineStart, 1.)));
   	l3p.x = mix(-1500., -200., expEasingIn(scene(t, lineStart, 1.)));
    enableL3 = t > lineStart;
    
    float rotateStart2 = lineStart + 1.1;
    l1Angle += mix(0., PI_2 , circEasingInOut(scene(t, rotateStart2, 1.)));
    l2Angle += mix(0., PI_2 , circEasingInOut(scene(t, rotateStart2, 1.)));
	l3Angle += mix(0., PI_2 , circEasingInOut(scene(t, rotateStart2, 1.)));    
	lpRotateAngle1 += mix(0., PI_2, circEasingInOut(scene(t, rotateStart2, 1.)));    
    l4p.y = mix(-1500., -200., expEasingIn(scene(t, rotateStart2, 1.)));
    enableL4 = t > rotateStart2;
    
    float rotateStart3 = rotateStart2 + 1.1;
    l1Angle += mix(0., PI_2 , circEasingInOut(scene(t, rotateStart3, 1.)));
    l2Angle += mix(0., PI_2 , circEasingInOut(scene(t, rotateStart3, 1.)));
	l3Angle += mix(0., PI_2 , circEasingInOut(scene(t, rotateStart3, 1.)));  
   	l4Angle += mix(0., PI_2 , circEasingInOut(scene(t, rotateStart3, 1.)));
    lpRotateAngle1 += mix(0., PI_2, circEasingInOut(scene(t, rotateStart3, 1.)));    
	float lpRotateAngle2 = mix(0., PI_2, circEasingInOut(scene(t, rotateStart3, 1.))); 

    float rotateStart4 = rotateStart3 + 1.3;
    l1Angle += mix(0., PI_2 , circEasingInOut(scene(t, rotateStart4, 1.)));
    l2Angle += mix(0., PI_2 , circEasingInOut(scene(t, rotateStart4, 1.)));
	l3Angle += mix(0., PI_2 , circEasingInOut(scene(t, rotateStart4, 1.)));  
   	l4Angle += mix(0., PI_2 , circEasingInOut(scene(t, rotateStart4, 1.)));
    lpRotateAngle1 += mix(0., PI_2, circEasingInOut(scene(t, rotateStart4, 1.)));    
	lpRotateAngle2 += mix(0., PI_2, circEasingInOut(scene(t, rotateStart4, 1.))); 
    scaleFactor += mix(0., -600., circEasingInOut(scene(t, rotateStart4, 1.)));

    float rotateStart5 = rotateStart4 + 1.1;
    l1Angle += mix(0., PI_2 , circEasingInOut(scene(t, rotateStart5, 1.)));
    l2Angle += mix(0., PI_2 , circEasingInOut(scene(t, rotateStart5, 1.)));
	l3Angle += mix(0., PI_2 , circEasingInOut(scene(t, rotateStart5, 1.)));  
   	l4Angle += mix(0., PI_2 , circEasingInOut(scene(t, rotateStart5, 1.)));
    lpRotateAngle1 += mix(0., PI_2, circEasingInOut(scene(t, rotateStart5, 1.)));    
	lpRotateAngle2 += mix(0., PI_2, circEasingInOut(scene(t, rotateStart5, 1.))); 
 
    float rotateStart6 = rotateStart5 + 1.3;
    l1Angle += mix(0., PI_2 , circEasingInOut(scene(t, rotateStart6, 1.)));
    l2Angle += mix(0., PI_2 , circEasingInOut(scene(t, rotateStart6, 1.)));
	l3Angle += mix(0., PI_2 , circEasingInOut(scene(t, rotateStart6, 1.)));  
   	l4Angle += mix(0., PI_2 , circEasingInOut(scene(t, rotateStart6, 1.)));
    lpRotateAngle1 += mix(0., PI_2, circEasingInOut(scene(t, rotateStart6, 1.)));    
	lpRotateAngle2 += mix(0., PI_2, circEasingInOut(scene(t, rotateStart6, 1.))); 
    
    float rotateStart7 = rotateStart6 + 1.3;
    l1Angle += mix(0., PI_2 , circEasingInOut(scene(t, rotateStart7, 1.)));
    l2Angle += mix(0., PI_2 , circEasingInOut(scene(t, rotateStart7, 1.)));
	l3Angle += mix(0., PI_2 , circEasingInOut(scene(t, rotateStart7, 1.)));  
   	l4Angle += mix(0., PI_2 , circEasingInOut(scene(t, rotateStart7, 1.)));
    lpRotateAngle1 += mix(0., PI_2, circEasingInOut(scene(t, rotateStart7, 1.)));    
	lpRotateAngle2 += mix(0., PI_2, circEasingInOut(scene(t, rotateStart7, 1.)));
    
    float rotateStart8 = rotateStart7 + 1.3;
    l1Angle += mix(0., PI_2 , circEasingInOut(scene(t, rotateStart8, 1.)));
    l2Angle += mix(0., PI_2 , circEasingInOut(scene(t, rotateStart8, 1.)));
	l3Angle += mix(0., PI_2 , circEasingInOut(scene(t, rotateStart8, 1.)));  
   	l4Angle += mix(0., PI_2 , circEasingInOut(scene(t, rotateStart8, 1.)));
    lpRotateAngle1 += mix(0., PI_2, circEasingInOut(scene(t, rotateStart8, 1.)));    
	lpRotateAngle2 += mix(0., PI_2, circEasingInOut(scene(t, rotateStart8, 1.)));

    float rotateStart9 = rotateStart8 + 1.3;
    l1Angle += mix(0., PI_2 , circEasingInOut(scene(t, rotateStart9, 1.)));
    l2Angle += mix(0., PI_2 , circEasingInOut(scene(t, rotateStart9, 1.)));
	l3Angle += mix(0., PI_2 , circEasingInOut(scene(t, rotateStart9, 1.)));  
   	l4Angle += mix(0., PI_2 , circEasingInOut(scene(t, rotateStart9, 1.)));
    lpRotateAngle1 += mix(0., PI_2, circEasingInOut(scene(t, rotateStart9, 1.)));    
	lpRotateAngle2 += mix(0., PI_2, circEasingInOut(scene(t, rotateStart9, 1.)));

    float rotateStart10 = rotateStart9 + 1.2;
    l1Angle += mix(0., PI_2 , circEasingInOut(scene(t, rotateStart10, 1.)));
    l2Angle += mix(0., PI_2 , circEasingInOut(scene(t, rotateStart10, 1.)));
	l3Angle += mix(0., PI_2 , circEasingInOut(scene(t, rotateStart10, 1.)));  
   	l4Angle += mix(0., PI_2 , circEasingInOut(scene(t, rotateStart10, 1.)));
    lpRotateAngle1 += mix(0., PI_2, circEasingInOut(scene(t, rotateStart10, 1.)));    
	lpRotateAngle2 += mix(0., PI_2, circEasingInOut(scene(t, rotateStart10, 1.)));
        
    
    float circleStart = rotateStart5;
    float deformR = 100.;
    c1.z += mix(0., deformR, expEasingOut(scene(t, rotateStart5, 1.)));
    c1.z += mix(0., -deformR, expEasingOut(scene(t, rotateStart6, 1.)));
	
	c3.z += mix(0., deformR, expEasingOut(scene(t, rotateStart6, 1.)));
  	c3.z += mix(0., -deformR, expEasingOut(scene(t, rotateStart7, 1.)));

    c5.z += mix(0., 50., circEasingInOut(scene(t, rotateStart7, 1.)));
    
    scaleFactor += mix(0., -300., circEasingInOut(scene(t, rotateStart7, 1.)));
    scaleFactor /= mix(1., 4., circEasingInOut(scene(t, rotateStart8, 1.)));
    scaleFactor /= mix(1., 4., circEasingInOut(scene(t, rotateStart9, 1.)));
    scaleFactor /= mix(1., 4., circEasingInOut(scene(t, rotateStart10, 1.)));

    float endingStart = rotateStart10 + 1.3;
    if(t >= endingStart){
        
        float rotateStart11 = endingStart;
   		l1Angle -= mix(0., PI_2 , circEasingInOut(scene(t, rotateStart11, 1.)));
    	l2Angle -= mix(0., PI_2 , circEasingInOut(scene(t, rotateStart11, 1.)));
		l3Angle -= mix(0., PI_2 , circEasingInOut(scene(t, rotateStart11, 1.)));  
   		l4Angle -= mix(0., PI_2 , circEasingInOut(scene(t, rotateStart11, 1.)));
    	lpRotateAngle1 -= mix(0., PI_2, circEasingInOut(scene(t, rotateStart11, 1.)));    
		lpRotateAngle2 -= mix(0., PI_2, circEasingInOut(scene(t, rotateStart11, 1.)));
        
        float rotateStart12 = rotateStart11 + 1.3;
   		l1Angle -= mix(0., PI_2 , circEasingInOut(scene(t, rotateStart12, 1.)));
    	l2Angle -= mix(0., PI_2 , circEasingInOut(scene(t, rotateStart12, 1.)));
		l3Angle -= mix(0., PI_2 , circEasingInOut(scene(t, rotateStart12, 1.)));  
   		l4Angle -= mix(0., PI_2 , circEasingInOut(scene(t, rotateStart12, 1.)));
    	lpRotateAngle1 -= mix(0., PI_2, circEasingInOut(scene(t, rotateStart12, 1.)));    
		lpRotateAngle2 -= mix(0., PI_2, circEasingInOut(scene(t, rotateStart12, 1.)));
        
        float rotateStart13 = rotateStart12+1.3;
   		l1Angle -= mix(0., PI_2 , circEasingInOut(scene(t, rotateStart13, 1.)));
    	l2Angle -= mix(0., PI_2 , circEasingInOut(scene(t, rotateStart13, 1.)));
		l3Angle -= mix(0., PI_2 , circEasingInOut(scene(t, rotateStart13, 1.)));  
   		l4Angle -= mix(0., PI_2 , circEasingInOut(scene(t, rotateStart13, 1.)));
    	lpRotateAngle1 -= mix(0., PI_2, circEasingInOut(scene(t, rotateStart13, 1.)));    
		lpRotateAngle2 -= mix(0., PI_2, circEasingInOut(scene(t, rotateStart13, 1.)));
        
        
        float endingTime = 1.;
   		
		c5.z  = mix(50., 0., circEasingInOut(scene(t, rotateStart12, endingTime)));
        l1p.x = mix(200., 5000., circEasingIn(scene(t, rotateStart13, endingTime)));  
		l2p.x = mix(-200., -5000., circEasingIn(scene(t, rotateStart13, endingTime)));
    	l3p.y = mix(0., -5000., circEasingIn(scene(t, rotateStart13, endingTime)));  
		l4p.y = mix(-200., -5000., circEasingIn(scene(t, rotateStart13, endingTime)));
        c1.z = mix(R, 0., circEasingInOut(scene(t, rotateStart13, endingTime)));
    	c2.z = mix(R, 0., circEasingInOut(scene(t, rotateStart13, endingTime)));
		c3.z = mix(R, 0., circEasingInOut(scene(t, rotateStart13, endingTime)));
		c4.z = mix(R, 0., circEasingInOut(scene(t, rotateStart13, endingTime)));

        theta += mix(0., -PI_2, circEasingInOut(scene(t, rotateStart13, 1.)));

        scaleFactor *= mix(1., 650., circEasingInOut(scene(t, endingStart, 1.)));

    }
    
    mat2 lpRotate1 = getRotationMat2(lpRotateAngle1);    
    l1m = getRotationMat2(l1Angle);
    l1mInv = getRotationMat2(-l1Angle);
    l2m = getRotationMat2(l2Angle);
    l2mInv = getRotationMat2(-l2Angle);
    l3m = getRotationMat2(l3Angle);
    l3mInv = getRotationMat2(-l3Angle);    
	l4m = getRotationMat2(l4Angle);
    l4mInv = getRotationMat2(-l4Angle); 
    
    l1p = lpRotate1 * l1p;
	l2p = lpRotate1 * l2p;
	l3p = lpRotate1 * l3p;
    
    mat2 lpRotate2 = getRotationMat2(lpRotateAngle2);    
 	l4p = lpRotate2 * l4p;    
    
    mat2 m = getRotationMat2(theta);    

    for(float i = 0. ; i < SAMPLE_NUM ; i++){
        vec2 position = ( (fragCoord.xy + rand2n(fragCoord.xy, i)) / iResolution.yy ) - vec2(ratio, 0.5);
        
        position = position * scaleFactor;
        position = m * position;
        
        
        float loopNum = IIS(position);

        if(loopNum >  0.){
            sum += hsv2rgb(vec3(0.01 + 0.05 * (loopNum-1.),1.0,1.0));
        }
    }
    fragColor = vec4(gammaCorrect(sum/SAMPLE_NUM), 1.);
}
