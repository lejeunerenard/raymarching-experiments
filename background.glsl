vec3 getBackground (in vec2 uv) {
  float coord = 1.0 * uv.y;

  return pow(#95E6FF, vec3(2.2)) * mix(0.2, 1.0, coord);
}
vec3 background = vec3(0.);
