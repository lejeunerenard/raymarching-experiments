#ifndef TWO_PI
#define TWO_PI 6.2831853072
#endif

#pragma glslify: rotMat2 = require(./rotation-matrix2)

vec3 dMin (vec3 d1, vec3 d2) {
  return (d1.x < d2.x) ? d1 : d2;
}

// Expects the following:
// - map : SDF w/ a vec3 for position & a float for the index
// - repetitions : SDF w/ a vec3 for position & a float for the index
vec3 pieSlice (in vec3 p) {
  vec3 d = vec3(1000000.0, 0, 0);

  float fRepetitions = float(repetitions);

  const float warpScale = 1.;

  float angle = 0.;
  float a = atan(p.y, p.x);

  for (int i = 0; i < repetitions; i++) {
    vec3 wQ = p;
    wQ.xy *= rotMat2(angle);
    vec3 b = map(wQ, float(i));
    d = dMin(d, b);
    angle += TWO_PI / fRepetitions;
  }

  return d;
}

#pragma glslify: export(pieSlice)
