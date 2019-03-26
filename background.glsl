vec3 getBackground (in vec2 uv) {
  // Convert from [-1,1] -> [0, 1]
  vec2 coord = 0.5 * (uv.xy + vec2(1.0));

  vec3 color = pow(#FF22C2, vec3(2.2));

  return color;
}
vec3 background = vec3(0.);
