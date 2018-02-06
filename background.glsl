vec3 getBackground (in vec2 uv) {
  // Convert from [-1,1] -> [0, 1]
  vec2 coord = 0.5 * (uv.xy + vec2(1.0));

  vec3 topColor = mix(#9EF4FF, #37CCB5, coord.x);
  vec3 color = mix(#FFBC5E, topColor, coord.y);

  color *= (0.97 + 0.05 * cnoise2(vec2(1.0 + cnoise2(uv), 33) * uv));

  return color;
}
vec3 background = vec3(0.);
