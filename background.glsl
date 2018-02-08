vec3 getBackground (in vec2 uv) {
  // Convert from [-1,1] -> [0, 1]
  vec2 coord = 0.5 * (uv.xy + vec2(1.0));
  float n = 0.2 * cnoise2(vec2(1.0, 13.0) * uv);
  n += 0.1 * cnoise2(vec2(0.5, 19.0) * uv + vec2(0.5, 1.0));


  vec3 color = mix(pow(#FF2E2E, vec3(2.2)), #FF7F12, 0.5 * n + coord.y);

  return color;
}
vec3 background = vec3(0.);
