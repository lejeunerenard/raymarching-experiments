#ifndef Iterations
  #define Iterations 10
#endif

float trapCalc (in vec3 p, in float k) {
  return dot(p, p) / (k * k);
}

#pragma glslify: octahedronFold = require(./folds/octahedron-fold, Iterations=Iterations, kifsM=kifsM, trapCalc=trapCalc)

vec2 kifs( inout vec3 p ) {
  float minD = 10000.;

  p = octahedronFold(p, minD);

  return vec2((length(p) - 1.) * pow(scale, - float(Iterations)), minD);
}

#pragma glslify: export(kifs)
