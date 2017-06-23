vec3 getBackground (in vec2 uv) {
  float coord = 1.0 * uv.y;

  // DARK Grey
  return pow(#9BCCFF, vec3(2.2));
}
vec3 background = vec3(0.);
