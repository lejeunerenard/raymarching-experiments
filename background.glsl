vec3 getBackground (in vec2 uv) {
  float coord = 1.0 * uv.y;

  return mix(#6284A3, #eeeefe, coord);
}
vec3 background = vec3(0.);
