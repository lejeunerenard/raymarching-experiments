#pragma glslify: pModPolarBack = require(./hg_sdf/p-mod-polar-c.glsl)
#pragma glslify: rotMat2Back = require(./rotation-matrix2)

vec3 getBackground (in vec2 uv) {
  // Convert from [-1,1] -> [0, 1]
  vec2 coord = 0.5 * (uv.xy + vec2(1.0));

  vec3 color = mix(vec3(0.7, 1, 0.6), vec3(1.0), saturate(length(coord)));
  color = mix(vec3(0.), vec3(0.05, 0.05, 0.1), saturate(length(coord)));
  color = mix(color, vec3(0.2, 0.2, 0.25), 0.2 * cnoise2(234. * uv) + smoothstep(0.1 + edge, 0.1, abs(uv.x)));

  return color;
}
vec3 background = vec3(0.);
