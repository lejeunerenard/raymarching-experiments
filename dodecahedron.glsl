#ifndef Iterations
  #define Iterations 10
#endif
const int IT = Iterations;

#pragma glslify: dodecahedronFold = require(./folds/dodecahedron-fold, Iterations=IT, kifsM=kifsM)

float sdBox( vec3 p, vec3 b ) {
  vec3 d = abs(p) - b;
  return min(max(d.x,max(d.y,d.z)),0.0) + length(max(d,0.0));
}

vec2 dodecasierpinski(in vec3 p) {
  float minD = 10000.;

  p = dodecahedronFold(p, minD);

  return vec2((sdBox(p, vec3(1.))) * pow(scale, -float(Iterations)), minD);
}

#pragma glslify: export(dodecasierpinski)
