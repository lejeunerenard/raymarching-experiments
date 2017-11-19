vec3 getBackground (in vec2 uv) {
  float coord = 1.0 * uv.y;

  return mix(vec3(0.8), vec3(1.00), coord);
}
vec3 background = vec3(0.);
