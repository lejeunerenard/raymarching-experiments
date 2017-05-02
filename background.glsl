vec3 getBackground (in vec2 uv) {
  return mix(#91E3B9, #80F8FF, 2.0 * uv.y - 1.0);
}
vec3 background = vec3(0.);
