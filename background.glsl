#pragma glslify: pModPolarBack = require(./hg_sdf/p-mod-polar-c.glsl)
#pragma glslify: rotMat2Back = require(./rotation-matrix2)
// #pragma glslify: snoise2 = require(glsl-noise/simplex/2d)

vec3 getBackground (in vec2 uv, in float universe) {
  // Convert from [-1,1] -> [0, 1]
  vec2 coord = 0.5 * (uv.xy + vec2(1.0));

  // vec3 color = mix(vec3(0.35), vec3(0.1), saturate(length(uv)));
  // vec3 color = mix(vec3(0.9), vec3(0.6), saturate(length(uv)));
  float i = saturate(pow(length(uv), 0.25));
  // i += snoise2(uv);
  vec3 color = mix(#FF3B93, vec3(0.05), i);

  return color;
}

vec3 getBackground (in vec2 uv) {
  return getBackground(uv, 0.);
}

vec3 background = vec3(0.);
