#ifndef Iterations
  #define Iterations 10
#endif

#pragma glslify: fold = require(./folds)
#pragma glslify: foldInv = require(./foldInv)

vec2 MengerSphere (inout vec3 p) {
  p = .5*p + vec3(.5);
  vec3 pp = abs(p -.5) - .5;

  float k = 2.0;

  float minD = 10000.;

  float d1 = max(pp.x, max(pp.y, pp.z));
  float d = d1;

  for (int i = 0; i < Iterations; i++) {
    vec3 pa = mod(3. * p * k, 3.);
    k *= scale;

    pp = .5 - abs(pa-1.5);

    pp = (vec4(pp, 1.) * kifsM).xyz;

    //distance inside the 3 axis aligned square tubes
    d1 = min(max(pp.x, pp.z), min(max(pp.x, pp.y), max(pp.y, pp.z))) / k;

    d = max(d, d1); // intersect

    minD = min(abs(dot(pp, pp)), minD);
  }

  float e = clamp(length(pp) - 10., 0., 100.);
  d = max(d, e);
  return vec2(d, minD);
}

#pragma glslify: export(MengerSphere)
