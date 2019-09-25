#pragma glslify: pModPolarBack = require(./hg_sdf/p-mod-polar-c.glsl)
#pragma glslify: rotMat2Back = require(./rotation-matrix2)

vec3 getBackground (in vec2 uv) {
  // Convert from [-1,1] -> [0, 1]
  vec2 coord = 0.5 * (uv.xy + vec2(1.0));

  const float size = 0.0075;

  vec3 color = mix(#668FCC, vec3(1), saturate(0.65 * length(coord)));
  color = #668FCC;

  return color;
}
vec3 background = vec3(0.);
