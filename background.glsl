vec3 getBackground (in vec2 uv) {
  float coord = 1.0 * uv.y;

  return pow(mix(#991307, #FF200C, coord), vec3(2.2));
}
vec3 background = vec3(0.);
