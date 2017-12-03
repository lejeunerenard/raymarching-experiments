vec3 getBackground (in vec2 uv) {
  vec2 coord = 2.0 * (uv.xy - vec2(0.5));

  return vec3(1);
}
vec3 background = vec3(0.);
