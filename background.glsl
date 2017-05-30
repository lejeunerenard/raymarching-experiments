vec3 getBackground (in vec2 uv) {
  return pow(mix(#C8CADB, #eeeeee, 0.5 * uv.y + 0.5), vec3(2.0));
}
vec3 background = vec3(0.);
