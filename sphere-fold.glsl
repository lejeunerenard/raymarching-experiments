void sphereFold (inout vec3 p, inout float dz) {
  const float fixedRadius2 = 1.;
  const float minRadius2 = .25;

  float r2 = dot(p, p);
  if (r2 < minRadius2) {
    float temp = (fixedRadius2 / minRadius2);
    p *= temp;
    dz *= temp;
  } else if (r2 < fixedRadius2) {
    float temp = (fixedRadius2 / r2);
    p *= temp;
    dz *= temp;
  }
} 

#pragma glslify: export(sphereFold)
