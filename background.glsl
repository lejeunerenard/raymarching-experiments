vec3 getBackground (in vec2 uv) {
  float coord = 1.0 * uv.y;

  return mix(#aaaaaa, #ffffff, coord);
}
vec3 background = vec3(0.);
