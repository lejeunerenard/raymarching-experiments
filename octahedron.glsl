#ifndef Iterations
  #define Iterations 10
#endif

float trapCalc (in vec3 p, in float k) {
  return dot(p, p) / (k * k);
}

float sdBox( vec3 p, vec3 b ) {
  vec3 d = abs(p) - b;
  return min(max(d.x,max(d.y,d.z)),0.0) + length(max(d,0.0));
}

#pragma glslify: octahedronFold = require(./folds/octahedron-fold, Iterations=Iterations, kifsM=kifsM, trapCalc=trapCalc)

vec2 kifs( inout vec3 p ) {
  float minD = 10000.;

  p = octahedronFold(p, minD);

  return vec2(sdBox(p, vec3(0.5)) * pow(scale, - float(Iterations)), minD);
}

#pragma glslify: export(kifs)
