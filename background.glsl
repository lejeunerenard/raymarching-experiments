#pragma glslify: pModPolarBack = require(./hg_sdf/p-mod-polar-c.glsl)
#pragma glslify: rotMat2Back = require(./rotation-matrix2)

vec3 getBackground (in vec2 uv) {
  // Convert from [-1,1] -> [0, 1]
  vec2 coord = 0.5 * (uv.xy + vec2(1.0));

  float l = length(uv);
  vec3 color = mix(vec3(0.8), vec3(1.0), l);
  color = vec3(0);
  // color = pow(color, vec3(0.25));

  // color = mix(color, vec3(0.75), 0.3);
  // color *= 0.8;

  return color;
}
vec3 background = vec3(0.);
