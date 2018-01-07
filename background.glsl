vec3 getBackground (in vec2 uv) {
  // Convert from [-1,1] -> [0, 1]
  vec2 coord = 0.5 * (uv.xy + vec2(1.0));

  const vec3 lightBlue = vec3(pow(#6328E8, vec3(2.2)));
  const vec3 superLightBlue = mix(lightBlue, vec3(0), 0.35);
  return mix(lightBlue, superLightBlue, coord.y);
}
vec3 background = vec3(0.);
