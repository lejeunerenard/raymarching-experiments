vec3 getBackground (in vec2 uv) {
  return mix(pow(#E6C6C5, vec3(2.2)), pow(#FFE6DF, vec3(2.2)), 0.5 * uv.y + 0.5);
}
vec3 background = vec3(0.);
