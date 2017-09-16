vec3 getBackground (in vec2 uv) {
  float coord = 1.0 * uv.y;

  return mix(#020000, #040202, coord);
}
vec3 background = vec3(0.);
