
#define POSITIVE_INF 1e+8
#define SKY_DIST     1e+7
#define STEP_MAX     1.0
#define MAX_STEPS    256
#define HIT_PREC     0.0001

float sceneDistFunction( float3 pos ){
    float dist = POSITIVE_INF;
    //===SCENE
    //printf( "%f %f %f %e \n", pos.x, pos.y, pos.z, dist );
    return dist;    
}

float castRay( float3 rd, float3 ro, float tmax ){
    float t = 0.0f;
    for( int i=0; i<MAX_STEPS; i++ ){
	    float dist = sceneDistFunction( ro+rd*t );
        if( dist<(HIT_PREC*t) ) break;
        if( t>tmax            ) return POSITIVE_INF;
        t += dist;
    }
    //printf("t %f == ro (%f,%f,%f) rd (%f,%f,%f) \n", t, ro.x,ro.y,ro.z, rd.x,rd.y,rd.z );
    return t;
}

float3 calcNormal( float3 pos ){
    /*
    float2 e = ((float2)(1.0f,-1.0f))*0.5773f*0.0005f;
    return normalize( 
        e.xyy*sceneDistFunction( pos + e.xyy ) + 
        e.yyx*sceneDistFunction( pos + e.yyx ) + 
        e.yxy*sceneDistFunction( pos + e.yxy ) + 
        e.xxx*sceneDistFunction( pos + e.xxx )
	);
	*/
	
	/*
	float2 e = (float2)(0.0f,0.001f);
	float v0 = sceneDistFunction( pos );
	return normalize( (float3)(
	    v0 - sceneDistFunction( pos+e.yxx ),
	    v0 - sceneDistFunction( pos+e.xyx ),
	    v0 - sceneDistFunction( pos+e.xxy )
	) );
	*/
	
	float2 e = (float2)(0.0f,0.001f);
	float3 nor = (float3)(
	    sceneDistFunction( pos+e.yxx ) - sceneDistFunction( pos-e.yxx ),
	    sceneDistFunction( pos+e.xyx ) - sceneDistFunction( pos-e.xyx ),
	    sceneDistFunction( pos+e.xxy ) - sceneDistFunction( pos-e.xxy )
	);
	return normalize(nor);
}

__kernel void rayTrace_basic(
    __global  float4* rayDir,    //  normalized_dir, step_length
    __global  float4* ray0,      //  pos, -
    __global  float4* hitPos,    //  pos, t(depth)
    __global  float4* hitNormal, //  normal, integral
    float max_depth
){
    
    float3 rd = rayDir[get_global_id(0)].xyz;
    float3 ro = ray0  [get_global_id(0)].xyz;
    //printf("== ro (%f,%f,%f) rd (%f,%f,%f) \n", ro.x,ro.y,ro.z, rd.x,rd.y,rd.z );
    //float t    = castRay( rd,ro, max_depth );
    
    float sum = 1.0f;
    float t   = 0.0f;
    for( int i=0; i<MAX_STEPS; i++ ){
	    float dist = sceneDistFunction( ro+rd*t );
	    //sum += 1.0f; 
	    //sum += 1.0f+dist; 
	    sum += 1.0f/(1.0f+80.0f*dist*dist); 
        if( dist<(HIT_PREC*t) ) break;
        if( t>max_depth       ){ t=POSITIVE_INF; break; };
        t += dist;
    }
    
    float3 pos = ro + t*rd;
    float3 nor = 0.0f;
    if( t<SKY_DIST ){
        nor = calcNormal( pos );  
        //printf( " == nor (%f,%f,%f) \n ", nor.x,nor.y,nor.z );
    };
    hitPos   [get_global_id(0)]     = (float4)(pos,t);
    hitNormal[get_global_id(0)]     = (float4)(nor,sum);
}


