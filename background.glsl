vec3 getBackground (in vec2 uv) {
  return mix(#FF67A6, #FF7114, 2.0 * uv.y - 1.0);
}
vec3 background = vec3(0.);
