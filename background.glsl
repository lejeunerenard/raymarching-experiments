#pragma glslify: pModPolarBack = require(./hg_sdf/p-mod-polar-c.glsl)
#pragma glslify: pMod2 = require(./hg_sdf/p-mod2.glsl)
#pragma glslify: rotMat2Back = require(./rotation-matrix2)
// #pragma glslify: snoise2 = require(glsl-noise/simplex/2d)

vec3 getBackground (in vec2 uv, in float universe) {
  // Convert from [-1,1] -> [0, 1]
  vec2 coord = 0.5 * (uv.xy + vec2(1.0));

  // float bgIndex = saturate(0.5 * (uv.y + 1.0));
  float bgIndex = saturate(length(uv));
  bgIndex = pow(bgIndex, 4.);
  // bgIndex = 1. - bgIndex;
  // vec3 color = mix(vec3(0.45, 0.3, 0.6), vec3(0.0), bgIndex);
  // vec3 color = mix(vec3(0.6), vec3(1.0), bgIndex);
  vec3 color = mix(vec3(1.0), #C8F7F3, bgIndex);
  // vec3 color = vec3(0);
  // color = mix(color, #FFC070, saturate(smoothstep(0.0, 0.5, uv.y)));

  return color;
}

vec3 getBackground (in vec2 uv) {
  return getBackground(uv, 0.);
}

vec3 background = vec3(0.);
