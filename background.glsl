vec3 getBackground (in vec2 uv) {
  // Convert from [-1,1] -> [0, 1]
  vec2 coord = 0.5 * (uv.xy + vec2(1.0));
  return vec3(0.);
  return mix(pow(#ABCAE8, vec3(2.2)), pow(#C9F2FF, vec3(2.2)), coord.y);
}
vec3 background = vec3(0.);
