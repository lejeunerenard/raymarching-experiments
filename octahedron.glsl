#ifndef Iterations
  #define Iterations 10
#endif

#pragma glslify: fold = require(./folds)
#pragma glslify: foldInv = require(./foldInv)

vec2 kifs( inout vec3 p ) {
  float minD = 10000.;

  for (int i = 0; i < Iterations; i++) {
    p = -abs(-p);

    // Folding
    fold(p.xy);
    fold(p.xz);
    foldInv(p.xy);
    foldInv(p.xz);

    // Stretch
    p = (vec4(p, 1.) * kifsM).xyz;

    minD = min(length(p), minD);
  }

  return vec2((length(p) - 1.) * pow(scale, - float(Iterations)), minD);
}

#pragma glslify: export(kifs)
