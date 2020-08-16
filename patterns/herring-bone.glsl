#ifndef PI
#define PI 3.1415926536
#define TWO_PI 6.2831853072
#endif

#pragma glslify: pMod1 = require(../hg_sdf/p-mod1.glsl)

float herringBone ( in vec2 q, float size, float borderThickness, float frequency, float stop ) {
  const float edge = 0.0025;

  float n = 1.;

  float c = pMod1(q.x, size);

  q.x = abs(q.x);

  n = dot(q, vec2(-1, 1));
  n *= frequency;
  n = sin(TWO_PI * n + phase(c));

  n = smoothstep(stop, stop + edge, n);
  n *= smoothstep(edge, 0., q.x - 0.5 * size + borderThickness);

  return n;
}

float herringBone ( in vec2 q, float size, float borderThickness, float frequency ) {
  return herringBone(q, size, borderThickness, frequency, 0.);
}

#pragma glslify: export(herringBone)
