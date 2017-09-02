vec3 getBackground (in vec2 uv) {
  float coord = 1.0 * uv.y;

  return mix(#334363, #889797, coord);
}
vec3 background = vec3(0.);
