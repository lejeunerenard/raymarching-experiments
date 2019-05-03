vec3 getBackground (in vec2 uv) {
  // Convert from [-1,1] -> [0, 1]
  vec2 coord = 0.5 * (uv.xy + vec2(1.0));

  vec3 colorOffset = vec3(0, 0.33, 0.73);
  vec2 absQ = abs(vec2(1.3, 1) * uv);
  float frameMask = smoothstep(edge, 0., max(absQ.x, absQ.y) - 0.5);
  vec3 color = 0.5 + 0.5 * cos(TWO_PI * (0.5 * coord.y + colorOffset));
  color *= smoothstep(0.5, 0., uv.y);
  vec3 other = vec3(0);
  color = mix(other, color, frameMask);
  return color;
}
vec3 background = vec3(0.);
