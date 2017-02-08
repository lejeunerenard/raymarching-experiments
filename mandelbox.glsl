void ballFold(inout vec3 v, inout float dz) {
  float r2 = dot(v, v);
  const float fixedRadius = 1.;
  const float fixedRadius2 = fixedRadius * fixedRadius;
  const float minRadius = 0.5;
  const float minRadius2 = minRadius * minRadius;

  if (r2 < minRadius2) {
    float temp = fixedRadius2 / minRadius2;
    v *= temp;
    dz *= temp;
  } else if ( r2 < fixedRadius2 ) {
    float temp = fixedRadius2 / r2;
    v *= temp;
    dz *= temp;
  }
}

#ifndef foldLimit
  #define foldLimit 1.
#endif

void boxFold (inout vec3 z) {
  z = clamp(z, -foldLimit, foldLimit) * 2. - z;
}

float minRadius2 = minRadius * minRadius;
const float fixedRadius = 0.9;
const float fixedRadius2 = fixedRadius * fixedRadius;
vec4 scalevec = vec4(s, s, s, abs(s)) / minRadius2;
float C1 = abs(s - 1.);
float C2 = pow(abs(s), float(1 - trap));

vec2 mandelbox (in vec3 z) {
  vec3 offset = z;
  float dr = 1.;
  vec4 p = vec4(z.xyz, 1.);
  vec4 p0 = vec4(z.xyz, 1.);

  float minD = maxDistance * 2.;

  for (int i = 0; i < trap; i++) {
    p.xyz = clamp(p.xyz, -1., 1.) * 2. - p.xyz;
    float r2 = dot(p.xyz, p.xyz);
    p.xyzw *= clamp(max(minRadius2/r2, minRadius2), 0., 1.);
    p *= rotM * max(mod(float(i), 4.), 1.) * rotM;

    p.xyzw = p*scalevec + p0;
    minD = min(minD, (length(p.xyz) - C1) / p.w);
  }

  return vec2((length(p.xyz) - C1) / p.w - C2, minD);
}

#pragma glslify: export(mandelbox)
