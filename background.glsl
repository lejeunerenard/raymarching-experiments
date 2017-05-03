vec3 getBackground (in vec2 uv) {
  return mix(#B839CC, #FFA0E8, 2.0 * uv.y - 1.0);
}
vec3 background = vec3(0.);
