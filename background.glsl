#pragma glslify: pModPolarBack = require(./hg_sdf/p-mod-polar-c.glsl)
#pragma glslify: rotMat2Back = require(./rotation-matrix2)

vec3 getBackground (in vec2 uv, in float universe) {
  // Convert from [-1,1] -> [0, 1]
  vec2 coord = 0.5 * (uv.xy + vec2(1.0));

  // vec3 color = mix(vec3(0.35), vec3(0.1), saturate(length(uv)));
  // vec3 color = mix(vec3(0.9), vec3(0.6), saturate(length(uv)));
  vec3 color = mix(#FF68AB, #E854E6, saturate(2. * uv.y));
  color *= 1.7;

  // color = mix(vec3(0.5), vec3(1.0), universe);

  return color;
}

vec3 getBackground (in vec2 uv) {
  return getBackground(uv, 0.);
}

vec3 background = vec3(0.);
