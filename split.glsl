#ifndef TWO_PI
#define TWO_PI 6.2831853072
#endif

// TODO properly test

float sdBox( vec3 p, vec3 b ) {
  vec3 d = abs(p) - b;
  return min(max(d.x,max(d.y,d.z)),0.0) + length(max(d,0.0));
}

#pragma glslify: pMod1 = require(./hg_sdf/p-mod1)

float split (inout vec3 p, in float splitTime) {
  float padding = size * (0.5 + 0.5 * cos(TWO_PI * splitTime - PI));
  float index = floor(p.y / (2.0 * padding + size));
  p.y = -index * 2.0 * padding + p.y;
  // Global offset
  p.y -= padding;

  vec3 cP = p;
  cP.y -= 0.5 * size + padding;
  float cI = pMod1(cP.y, 2.0 * padding + size);

  #ifndef breath
  const float breadth = 10.0;
  #endif

  return sdBox(cP, vec3(breadth, 0.5 * size, breadth));
}

#pragma glslify: export(split)
