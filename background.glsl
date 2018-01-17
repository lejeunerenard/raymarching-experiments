vec3 getBackground (in vec2 uv) {
  // Convert from [-1,1] -> [0, 1]
  vec2 coord = 0.5 * (uv.xy + vec2(1.0));

  vec3 defaultColor = mix(vec3(0.7), vec3(0.8), coord.y);

  vec3 color = mix(defaultColor, #33AAFF, smoothstep(0.95, 0., length(uv)));

  return color;
}
vec3 background = vec3(0.);
