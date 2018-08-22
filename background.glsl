vec3 getBackground (in vec2 uv) {
  // Convert from [-1,1] -> [0, 1]
  vec2 coord = 0.5 * (uv.xy + vec2(1.0));
  float n = -coord.y + 0.05 * cnoise2(345. * uv);
  n = pow(n, 1.5);
  return mix(pow(#D93741, vec3(2.2)), pow(#040404, vec3(2.2)), 1. - n);
}
vec3 background = vec3(0.);
