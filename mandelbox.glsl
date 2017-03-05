#ifndef foldLimit
  #define foldLimit 1.
#endif

#pragma glslify: fold = require(./folds)

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

  float k = 1.0;

  for (int i = 0; i < trap; i++) {
    p.xyz = clamp(p.xyz, -foldLimit, foldLimit) * 2. - p.xyz;

    float r2 = dot(p.xyz, p.xyz);
    p.xyzw *= clamp(max(minRadius2/r2, minRadius2), 0., 1.);

    p *= rotM * max(mod(float(i), 4.), 1.) * rotM;

    p.xyzw = p*scalevec + p0;
    minD = min(minD, (dot(p.xyz, p.xyz) - C1) / (k * k * p.w));
    k *= 1.1;
  }

  return vec2((length(p.xyz) - C1) / p.w - C2, minD);
}

#pragma glslify: export(mandelbox)
