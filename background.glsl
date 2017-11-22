vec3 getBackground (in vec2 uv) {
  float coord = 1.0 * uv.y;

  return mix(vec3(0.7), vec3(1), coord);
}
vec3 background = vec3(0.);
