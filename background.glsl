#pragma glslify: pModPolarBack = require(./hg_sdf/p-mod-polar-c.glsl)
#pragma glslify: pMod2 = require(./hg_sdf/p-mod2.glsl)
#pragma glslify: rotMat2Back = require(./rotation-matrix2)
// #pragma glslify: snoise2 = require(glsl-noise/simplex/2d)

vec3 getBackground (in vec2 uv, in float universe) {
  // Convert from [-1,1] -> [0, 1]
  vec2 coord = 0.5 * (uv.xy + vec2(1.0));

  // vec3 color = mix(#1C0036, vec3(0.0), saturate(length(uv)));
  vec3 color = mix(vec3(0.8), vec3(1.0), saturate(length(uv)));
  // vec3 color = mix(#815DFF, #FA69D8, saturate(uv.y + 0.5));

  return color;
}

vec3 getBackground (in vec2 uv) {
  return getBackground(uv, 0.);
}

vec3 background = vec3(0.);
