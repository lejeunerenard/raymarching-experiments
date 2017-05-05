vec3 getBackground (in vec2 uv) {
  return mix(#CCFCF3, #F0FEFF, 2.0 * uv.y - 1.0);
}
vec3 background = vec3(0.);
