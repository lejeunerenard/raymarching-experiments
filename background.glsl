#pragma glslify: pModPolar = require(./hg_sdf/p-mod-polar-c.glsl)
const float backgroundR = 0.5;
float backgroundMask (in vec2 uv, in float r) {
  pModPolar(uv, 6.);
  return smoothstep(edge, 0., uv.x - r);
}

vec3 getBackground (in vec2 uv) {
  // Convert from [-1,1] -> [0, 1]
  vec2 coord = 0.5 * (uv.xy + vec2(1.0));

  float frameMask = backgroundMask(uv, backgroundR);
  vec2 backgroundUv = uv;
  pModPolar(backgroundUv, 6.);
  vec3 colorOffset = vec3(0, 0.33, 0.73);
  vec3 color = vec3(smoothstep(0.5 * edge, 0., abs(backgroundUv.x - (backgroundR + 2. * edge))));
  return color;
}
vec3 background = vec3(0.);
