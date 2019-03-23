vec3 getBackground (in vec2 uv) {
  // Convert from [-1,1] -> [0, 1]
  vec2 coord = 0.5 * (uv.xy + vec2(1.0));

  const vec3 blue = pow(#3B4BFF, vec3(2.2));
  vec3 color = mix(pow(#1724B2, vec3(2.2)), blue, coord.y);
  color = blue;

  return color;
}
vec3 background = vec3(0.);
