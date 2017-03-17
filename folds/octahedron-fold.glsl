#ifndef Iterations
  #define Iterations 10
#endif

#pragma glslify: fold = require(../folds)
#pragma glslify: foldInv = require(../foldInv)

vec3 octahedronFold (in vec3 p, inout float minD) {
  float k = 1.0;

  for (int i = 0; i < Iterations; i++) {
    p = -abs(-p);

    // Folding
    fold(p.xy);
    fold(p.xz);
    foldInv(p.xy);
    foldInv(p.xz);

    // Stretch
    p = (vec4(p, 1.) * kifsM).xyz;

    minD = min(trapCalc(p, k), minD);
    k *= 1.1;
  }

  return p;
}

#pragma glslify: export(octahedronFold)
