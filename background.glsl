#pragma glslify: pModPolarBack = require(./hg_sdf/p-mod-polar-c.glsl)
#pragma glslify: rotMat2Back = require(./rotation-matrix2)

vec3 getBackground (in vec2 uv) {
  // Convert from [-1,1] -> [0, 1]
  vec2 coord = 0.5 * (uv.xy + vec2(1.0));

  vec3 color = vec3(0);
  // vec3 color = mix(vec3(0.7), vec3(1.0), saturate(length(uv)));
  // vec3 color = mix(vec3(0.0), vec3(0.2), saturate(length(uv)));

  return color;
}
vec3 background = vec3(0.);
