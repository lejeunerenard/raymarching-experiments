#ifndef Iterations
  #define Iterations 10
#endif

float trapCalc (in vec3 p, in float k) {
  return dot(p, p) / (k * k);
}

#pragma glslify: dodecahedronFold = require(./folds/dodecahedron-fold, Iterations=1, kifsM=kifsM)
#pragma glslify: octahedronFold = require(./folds/octahedron-fold, Iterations=1, kifsM=kifsM, trapCalc=trapCalc)

vec2 dodecasierpinski(in vec3 p) {
  float minD = 10000.;

  for (int i = 0; i < Iterations; i++) {
    p = dodecahedronFold(p, minD);
    p = octahedronFold(p, minD);
  }

  return vec2((sqrt(length(p)) - 2.0) * pow(scale, -float(Iterations * 2)), minD);
}

#pragma glslify: export(dodecasierpinski)
