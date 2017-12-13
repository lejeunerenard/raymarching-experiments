vec3 getBackground (in vec2 uv) {
  vec2 coord = 2.0 * (uv.xy - vec2(0.5));

  return mix(vec3(0.720, 0.800, 0.800), vec3(1.00, 1.00, 1.00), coord.y);
}
vec3 background = vec3(0.);
