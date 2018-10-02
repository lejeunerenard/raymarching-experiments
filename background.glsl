vec3 getBackground (in vec2 uv) {
  // Convert from [-1,1] -> [0, 1]
  vec2 coord = 0.5 * (uv.xy + vec2(1.0));

  // return vec3(max(abs(uv.x), abs(uv.y)));
  vec3 color = mix(vec3(0.), vec3(0.01), smoothstep(0., 1., length(uv)));
  color *= 0.5 + 0.5 * coord.y;
  return color;
}
vec3 background = vec3(0.);
