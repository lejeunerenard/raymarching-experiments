#pragma glslify: pModPolarBack = require(./hg_sdf/p-mod-polar-c.glsl)
#pragma glslify: rotMat2Back = require(./rotation-matrix2)

vec3 getBackground (in vec2 uv) {
  // Convert from [-1,1] -> [0, 1]
  vec2 coord = 0.5 * (uv.xy + vec2(1.0));

  // vec3 color = mix(vec3(0.35), vec3(0.0), saturate(length(uv)));
  vec3 color = mix(vec3(0.9), vec3(1.0), saturate(length(uv)));
  vec3 dI = vec3(dot(uv, 0.6 * vec2(0.2, 1.)) + 0.4);
  color = mix(#FFA742, #FF5074, smoothstep(0., 0.33, dI.x));
  color = mix(color, #D43CE8, smoothstep(0.33, 0.67, dI.x));
  color = mix(color, #8142FF, smoothstep(0.67, 1.00, dI.x));
  // color = vec3(0);

  return color;
}
vec3 background = vec3(0.);
