vec3 getBackground (in vec2 uv) {
  return mix(#4EFF7E, #CBFFAD, 2.0 * uv.y - 1.0);
}
vec3 background = vec3(0.);
