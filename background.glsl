vec3 getBackground (in vec2 uv) {
  float coord = 1.0 * uv.y;

  return mix(#666666, #aaaaaa, coord);
}
vec3 background = vec3(0.);
