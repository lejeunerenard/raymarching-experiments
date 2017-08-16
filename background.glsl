vec3 getBackground (in vec2 uv) {
  float coord = 1.0 * uv.y;

  return mix(#FF50CC, #CC40A3, coord);
}
vec3 background = vec3(0.);
