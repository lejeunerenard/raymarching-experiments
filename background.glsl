#pragma glslify: pModPolarBack = require(./hg_sdf/p-mod-polar-c.glsl)
#pragma glslify: rotMat2Back = require(./rotation-matrix2)

vec3 getBackground (in vec2 uv) {
  // Convert from [-1,1] -> [0, 1]
  vec2 coord = 0.5 * (uv.xy + vec2(1.0));

  const float size = 0.0075;
  vec2 c = floor(coord / size) / size;
  vec3 color = mix(vec3(0.), vec3(0.2), mod(dot(c, vec2(1)), 2.));
  color = vec3(0.1, 0.5785234899328859, 0.9);

  return color;
}
vec3 background = vec3(0.);
