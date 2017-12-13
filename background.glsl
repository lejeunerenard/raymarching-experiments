vec3 getBackground (in vec2 uv) {
  vec2 coord = 2.0 * (uv.xy - vec2(0.5));

  return mix(vec3(0.874), vec3(1), coord.y);
}
vec3 background = vec3(0.);
