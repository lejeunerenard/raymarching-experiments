vec3 getNormal( in vec3 p, in float eps ) {
    vec2 e = vec2(1.0,-1.0)*.05773*eps;
    return normalize(
      e.xyy * map( p + e.xyy ).x + 
      e.yyx * map( p + e.yyx ).x + 
      e.yxy * map( p + e.yxy ).x + 
      e.xxx * map( p + e.xxx ).x );
}

#pragma glslify: export(getNormal)
