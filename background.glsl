vec3 getBackground (in vec2 uv) {
  return mix(#000001, #020202, 0.5 * uv.y + 0.5);
}
vec3 background = vec3(0.);
