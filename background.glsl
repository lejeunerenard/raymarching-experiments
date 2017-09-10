vec3 getBackground (in vec2 uv) {
  float coord = 1.0 * uv.y;

  return mix(#aaaaaa, #cccccc, coord);
}
vec3 background = vec3(0.);
