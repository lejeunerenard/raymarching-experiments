vec3 getBackground (in vec2 uv) {
  float coord = 1.0 * uv.y;

  return mix(#8888a8, #ccccdc, coord);
}
vec3 background = vec3(0.);
