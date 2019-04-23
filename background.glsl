vec3 getBackground (in vec2 uv) {
  // Convert from [-1,1] -> [0, 1]
  vec2 coord = 0.5 * (uv.xy + vec2(1.0));

  float n = dot(coord, 57. * vec2(1));
  n = sin(n);
  n = smoothstep(edge, 0., n);

  vec3 color = vec3(n);
  return color;
}
vec3 background = vec3(0.);
