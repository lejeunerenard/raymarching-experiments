vec3 getBackground (in vec2 uv) {
  // Convert from [-1,1] -> [0, 1]
  vec2 coord = 0.5 * (uv.xy + vec2(1.0));

  vec3 color = mix(pow(#FF9D5C, vec3(2.2)), pow(#FF64D8, vec3(2.2)), coord.y);
  // color = pow(#C2EAFF, vec3(2.2));
  return color;
}
vec3 background = vec3(0.);
