vec3 getBackground (in vec2 uv) {
  return mix(#5a4444, #d9d9d9, pow(0.5 * uv.y + 0.5, 0.5));
}
vec3 background = vec3(0.);
