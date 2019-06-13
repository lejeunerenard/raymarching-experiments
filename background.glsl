#pragma glslify: pModPolar = require(./hg_sdf/p-mod-polar-c.glsl)
#pragma glslify: rotMat2 = require(./rotation-matrix2)
const float backgroundR = 0.5;
float backgroundMask (in vec2 uv, in float r) {
  return smoothstep(edge, 0., length(uv) - r);
}

vec3 getBackground (in vec2 uv) {
  // Convert from [-1,1] -> [0, 1]
  vec2 coord = 0.5 * (uv.xy + vec2(1.0));

  float frameMask = backgroundMask(uv, backgroundR);
  vec2 backgroundUv = uv;
  float l = length(uv);
  vec3 color = vec3(0.60 + 0.60 * saturate(1. - l));
  return color;
}
vec3 background = vec3(0.);
