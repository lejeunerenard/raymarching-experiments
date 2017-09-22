vec3 getBackground (in vec2 uv) {
  float coord = 1.0 * uv.y;

  return mix(#C4B3FF, #989CE8, coord);
}
vec3 background = vec3(0.);
