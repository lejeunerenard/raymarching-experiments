vec3 getBackground (in vec2 uv) {
  // Convert from [-1,1] -> [0, 1]
  vec2 coord = 0.5 * (uv.xy + vec2(1.0));

  const vec3 orange = #F7C09E;
  const vec3 beige = #FFD6BD;
  const vec3 darkGrey = #578071;
  const vec3 teal = #62FFC7;
  vec3 color = beige * 0.9;
  color = vec3(1);
  return color;
}
vec3 background = vec3(0.);
