vec3 getBackground (in vec2 uv) {
  // Convert from [-1,1] -> [0, 1]
  vec2 coord = 0.5 * (uv.xy + vec2(1.0));

  vec3 color = mix(vec3(0.45), vec3(0.80), length(coord));
  color = #5A9999;
  // color = #12CACC;
  color = #B8E0E0;

  return color;
}
vec3 background = vec3(0.);
