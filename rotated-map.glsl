#ifndef PI
#define PI 3.1415926536
#define TWO_PI 6.2831853072
#endif

#pragma glslify: rotationMatrix = require(./rotation-matrix3)

float rotatedMap (in vec3 p) {
  const float angleOffset = (TWO_PI + PI * 0.125) * 2.0 / float(rotSymmetry);
  mat3 rotOffset = rotationMatrix(vec3(0.0, 0.0, 1.0), angleOffset);
  mat3 modelRotOffset = rotationMatrix(vec3(0.0, 0.0, 1.0), 0.2 * TWO_PI);

  float outD = 1000.0;

  for (int i = 0; i < rotSymmetry; i++) {
    vec3 r = p;
    r.z += 0.25 * sin(PI * (0.75 * float(i)));

    float d = map(r, i);
    outD = min(outD, d);

    p *= rotOffset;
    p.xy *= 0.99;
  }

  return outD;
}

#pragma glslify: export(rotatedMap)
