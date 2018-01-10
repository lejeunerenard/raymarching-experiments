vec3 getBackground (in vec2 uv) {
  // Convert from [-1,1] -> [0, 1]
  vec2 coord = 0.5 * (uv.xy + vec2(1.0));

  const vec3 orange = pow(#BCACB3, vec3(2.2));
  return mix(orange, vec3(1), 0.8);
}
vec3 background = vec3(0.);
