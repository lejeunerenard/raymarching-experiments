vec3 getBackground (in vec2 uv) {
  // Convert from [-1,1] -> [0, 1]
  vec2 coord = 0.5 * (uv.xy + vec2(1.0));

  return (0.875 + 0.125 * sin(PI * 0.2 * time)) * mix(vec3(pow(#FFC399, vec3(2.2))), pow(#FFFFA9, vec3(2.2)), 0.75 * coord.y);
}
vec3 background = vec3(0.);
