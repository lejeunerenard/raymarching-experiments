vec3 getBackground (in vec2 uv) {
  // Convert from [-1,1] -> [0, 1]
  vec2 coord = 0.5 * (uv.xy + vec2(1.0));

  return mix(#2a66ff, #4a86ff, coord.y);
}
vec3 background = vec3(0.);
