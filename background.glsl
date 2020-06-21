#pragma glslify: pModPolarBack = require(./hg_sdf/p-mod-polar-c.glsl)
#pragma glslify: rotMat2Back = require(./rotation-matrix2)

vec3 getBackground (in vec2 uv, in float universe) {
  // Convert from [-1,1] -> [0, 1]
  vec2 coord = 0.5 * (uv.xy + vec2(1.0));

  // vec3 color = mix(vec3(0.35), vec3(0.1), saturate(length(uv)));
  // vec3 color = mix(vec3(0.9), vec3(0.6), saturate(length(uv)));
  // vec3 color = mix(#FFADC3, #E0E7FF, saturate(dot(uv, 1.5 * vec2(0.2, 0.8)) + 0.2));
  vec3 color = mix(#DCABFF, #ABE8FF, saturate(length(uv)));
  color *= 0.8;

  // vec3 color = vec3(0.1);

  // color = mix(vec3(0.5), vec3(1.0), universe);

  return color;
}

vec3 getBackground (in vec2 uv) {
  return getBackground(uv, 0.);
}

vec3 background = vec3(0.);
