bool sphere(vec3 org, vec3 dir, out float near, out float far)
{
	float b = dot(dir, org);
	float c = dot(org, org) - 0.25;
	float delta = b*b - c;
	if( delta < 0.0) 
		return false;
	float deltasqrt = sqrt(delta);
	near = -b - deltasqrt;
	far = -b + deltasqrt;
	return far > 0.0;
}

bool cylinder(vec3 org, vec3 dir, out float near, out float far)
{
	// quadratic x^2 + y^2 = 0.5^2 => (org.x + t*dir.x)^2 + (org.y + t*dir.y)^2 = 0.5
	float a = dot(dir.xy, dir.xy);
	float b = dot(org.xy, dir.xy);
	float c = dot(org.xy, org.xy) - 0.25;

	float delta = b * b - a * c;
	if( delta < 0.0 )
		return false;

	// 2 roots
	float deltasqrt = sqrt(delta);
	float arcp = 1.0 / a;
	near = (-b - deltasqrt) * arcp;
	far = (-b + deltasqrt) * arcp;
	
	// order roots
	float temp = min(far, near);
	far = max(far, near);
	near = temp;

	float znear = org.z + near * dir.z;
	float zfar = org.z + far * dir.z;

	// top, bottom
	vec2 zcap = vec2(0.5, -0.5);
	vec2 cap = (zcap - org.z) / dir.z;

	if ( znear < zcap.y )
		near = cap.y;
	else if ( znear > zcap.x )
		near = cap.x;

	if ( zfar < zcap.y )
		far = cap.y;
	else if ( zfar > zcap.x )
		far = cap.x;
	
	return far > 0.0 && far > near;
}

// cone inscribed in a unit cube centered at 0
bool cone(vec3 org, vec3 dir, out float near, out float far)
{
	// scale and offset into a unit cube
	org.x += 0.5;
	float s = 0.5;
	org.x *= s;
	dir.x *= s;
	
	// quadratic x^2 = y^2 + z^2
	float a = dir.y * dir.y + dir.z * dir.z - dir.x * dir.x;
	float b = org.y * dir.y + org.z * dir.z - org.x * dir.x;
	float c = org.y * org.y + org.z * org.z - org.x * org.x;
	
	float cap = (s - org.x) / dir.x;
	
	// linear
	if( a == 0.0 )
	{
		near = -0.5 * c/b;
		float x = org.x + near * dir.x;
		if( x < 0.0 || x > s )
			return false; 

		far = cap;
		float temp = min(far, near); 
		far = max(far, near);
		near = temp;
		return far > 0.0;
	}

	float delta = b * b - a * c;
	if( delta < 0.0 )
		return false;

	// 2 roots
	float deltasqrt = sqrt(delta);
	float arcp = 1.0 / a;
	near = (-b - deltasqrt) * arcp;
	far = (-b + deltasqrt) * arcp;
	
	// order roots
	float temp = min(far, near);
	far = max(far, near);
	near = temp;

	float xnear = org.x + near * dir.x;
	float xfar = org.x + far * dir.x;

	if( xnear < 0.0 )
	{
		if( xfar < 0.0 || xfar > s )
			return false;
		
		near = far;
		far = cap;
	}
	else if( xnear > s )
	{
		if( xfar < 0.0 || xfar > s )
			return false;
		
		near = cap;
	}
	else if( xfar < 0.0 )
	{
		// The apex is problematic,
		// additional checks needed to
		// get rid of the blinking tip here.
		far = near;
		near = cap;
	}
	else if( xfar > s )
	{
		far = cap;
	}
	
	return far > 0.0;
}

// cube() by Simon Green
bool cube(vec3 org, vec3 dir, out float near, out float far)
{
	// compute intersection of ray with all six bbox planes
	vec3 invR = 1.0/dir;
	vec3 tbot = invR * (-0.5 - org);
	vec3 ttop = invR * (0.5 - org);
	
	// re-order intersections to find smallest and largest on each axis
	vec3 tmin = min (ttop, tbot);
	vec3 tmax = max (ttop, tbot);
	
	// find the largest tmin and the smallest tmax
	vec2 t0 = max(tmin.xx, tmin.yz);
	near = max(t0.x, t0.y);
	t0 = min(tmax.xx, tmax.yz);
	far = min(t0.x, t0.y);

	// check for hit
	return near < far && far > 0.0;
}

// frustum inscribed in a unit cube centered at 0
#define INF 1.0e38
bool frustum(vec3 org, vec3 dir, float apex, out float near, out float far)
{
	vec2 dirf = vec2(0.5 - apex, 0.5); 
	vec3 tbot, ttop;
	
	// intersection with near and far planes
	float invdirx = 1.0 / dir.x;
	tbot.x = invdirx * (-0.5 - org.x);
	ttop.x = invdirx * (0.5 - org.x);

	float temp = dirf.y * (org.x-apex);
	
	// intersection with inclined planes on y
	tbot.y = (-temp - dirf.x * org.y) / (dirf.x * dir.y + dirf.y * dir.x);
	ttop.y = ( temp - dirf.x * org.y) / (dirf.x * dir.y - dirf.y * dir.x);
	
	// intersection with inclined planes on z
	tbot.z = (-temp - dirf.x * org.z) / (dirf.x * dir.z + dirf.y * dir.x);
	ttop.z = ( temp - dirf.x * org.z) / (dirf.x * dir.z - dirf.y * dir.x);
	
	// if intersecting behind the apex, set t to ray's end
	vec4 tempt = vec4(tbot.yz, ttop.yz);
	tempt = mix(tempt, INF * sign(dir.xxxx), step(org.xxxx + tempt * dir.xxxx, vec4(apex)));
	tbot.yz = tempt.xy;
	ttop.yz = tempt.zw;

	// re-order intersections to find smallest and largest on each axis
	vec3 tmin = min(ttop, tbot);
	vec3 tmax = max(ttop, tbot);
	
	// find the largest tmin and the smallest tmax
	vec2 t0 = max(tmin.xx, tmin.yz);
	near = max(t0.x, t0.y);
	t0 = min(tmax.xx, tmax.yz);
	far = min(t0.x, t0.y);

	// check for hit
	return near < far && far > 0.0;
}

void transformray (vec3 ro, vec3 rd, mat2 rotationY, vec3 offset, out vec3 outro, out vec3 outrd)
{
	outro = ro + offset;
	outro = vec3(rotationY * outro.xz, outro.y).xzy;
	outrd = vec3(rotationY * rd.xz, rd.y).xzy;
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
	// camera
	vec2 q = fragCoord.xy/iResolution.xy;
	vec2 p = -1.0 + 2.0 * q;
	p.x *= iResolution.x/iResolution.y;
	vec3 camro = normalize(vec3(1.0, -0.1, 0.0));
	vec3 w = -camro;
	camro *= 2.5;
	vec3 u = normalize(cross( vec3(0.0, 1.0, 0.0), w ));
	vec3 v = normalize(cross(w,u));
	vec3 camrd = normalize(p.x * u + p.y * v + 1.5 * w);
	fragColor = vec4(0.0);
	
	// rotation
	float angle = 5.0 * iMouse.x / iResolution.x;
	if( iMouse.z < 0.5 )
		angle = iTime + 4.7;
	float ca = cos(angle);
	float sa = sin(angle);
	mat2  m = mat2(ca, -sa, sa, ca);
	
	float far, near;
	vec3 ro, rd;
	
	// cube
	transformray(camro, camrd, m, vec3(0, -0.7, 0.8), ro, rd );
	if (cube (ro, rd, near, far))
		fragColor += vec4(far - max(near, 0.0));
	
	// frustum
	transformray(camro, camrd, m, vec3(0, -0.7, -0.8), ro, rd);
	if (frustum (ro, rd, -1.0, near, far))
		fragColor += vec4(far - max(near, 0.0));
	
	// sphere
	transformray(camro, camrd, m, vec3(0, 0.7, 1.5), ro, rd );
	if (sphere (ro, rd, near, far))
		fragColor += vec4(far - max(near, 0.0));

	// cylinder
	transformray(camro, camrd, m, vec3(0, 0.7, 0.0), ro, rd );
	if (cylinder (ro, rd, near, far))
		fragColor += vec4(far - max(near, 0.0));

	// cone
	transformray(camro, camrd, m, vec3(0, 0.7, -1.5), ro, rd);
	if (cone (ro, rd, near, far))
		fragColor += vec4(far - max(near, 0.0));
}