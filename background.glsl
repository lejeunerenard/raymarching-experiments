vec3 getBackground (in vec2 uv) {
  // Convert from [-1,1] -> [0, 1]
  vec2 coord = 0.5 * (uv.xy + vec2(1.0));

  vec3 color = mix(vec3(0.4), vec3(0.9), coord.y);
  const vec3 blue = pow(#1A11FF, vec3(2.2));
  const vec3 cyan = pow(#34EEFF, vec3(2.2));
  color = cyan;

  return color;
}
vec3 background = vec3(0.);
