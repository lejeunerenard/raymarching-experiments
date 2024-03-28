vec3 getNormal( in vec3 p, in float eps, in float generalT ) {
    vec2 e = vec2(1.0,-1.0)*.05773*eps;
    return normalize(
      e.xyy * map( p + e.xyy, generalT).x + 
      e.yyx * map( p + e.yyx, generalT).x + 
      e.yxy * map( p + e.yxy, generalT).x + 
      e.xxx * map( p + e.xxx, generalT).x );
}

#pragma glslify: export(getNormal)
