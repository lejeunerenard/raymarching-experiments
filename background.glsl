#pragma glslify: pModPolarBack = require(./hg_sdf/p-mod-polar-c.glsl)
#pragma glslify: rotMat2Back = require(./rotation-matrix2)

vec3 getBackground (in vec2 uv) {
  // Convert from [-1,1] -> [0, 1]
  vec2 coord = 0.5 * (uv.xy + vec2(1.0));

  float l = length(uv);
  vec3 color = mix(1.3 * #FFB0DD, vec3(1.0), saturate(uv.y));

  // color *= 0.8;

  return color;
}
vec3 background = vec3(0.);
