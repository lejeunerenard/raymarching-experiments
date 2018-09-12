vec3 getBackground (in vec2 uv) {
  // Convert from [-1,1] -> [0, 1]
  vec2 coord = 0.5 * (uv.xy + vec2(1.0));

  vec3 color = vec3(1);

  vec2 absQ = abs(uv + vec2(0.4, -0.45));
  color = mix(color, vec3(0), smoothstep(0.01, 0.0, max(absQ.x, absQ.y) - 0.35));

  absQ = abs(uv + vec2(0.1, 0.2));
  color = mix(color, vec3(1, 0, 0), smoothstep(0.01, 0.0, max(absQ.x, absQ.y) - 0.4));

  absQ = abs(uv - vec2(0.45, 0.35));
  color = mix(color, vec3(0, 0, 1), smoothstep(0.01, 0.0, max(absQ.x, absQ.y) - 0.3));

  return color;
}
vec3 background = vec3(0.);
