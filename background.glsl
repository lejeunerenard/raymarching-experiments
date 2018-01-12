vec3 getBackground (in vec2 uv) {
  // Convert from [-1,1] -> [0, 1]
  vec2 coord = 0.5 * (uv.xy + vec2(1.0));

  vec2 aUV = abs(uv);
  float i = smoothstep(0.0, 0.75, length(aUV));
  const vec3 purple = pow(#9636FF, vec3(2.2));
  return mix(purple, vec3(0.001), pow(i, 0.1));
}
vec3 background = vec3(0.);
