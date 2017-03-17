#ifndef foldLimit
  #define foldLimit 1.
#endif

float C1 = abs(s - 1.);
float C2 = pow(abs(s), float(1 - trap));

float trapCalc (in vec4 z) {
  return (length(z.xyz) - C1) / z.w;
}

#pragma glslify: fold = require(./folds)
#pragma glslify: mandelboxFold = require(./folds/mandelbox-fold, foldLimit=foldLimit, trap=trap, C1=C1, minRadius=minRadius, s=s, rotM=rotM, trapCalc=trapCalc)

vec2 mandelbox (in vec3 z) {
  vec4 p = vec4(z.xyz, 1.);

  float minD = maxDistance * 2.;

  p = mandelboxFold(p, minD);

  return vec2((length(p.xyz) - C1) / p.w - C2, minD);
}

#pragma glslify: export(mandelbox)
