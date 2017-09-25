vec3 getBackground (in vec2 uv) {
  float coord = 1.0 * uv.y;

  return mix(pow(#0AA4EA, vec3(2.2)), pow(#F87717, vec3(2.2)), coord);
}
vec3 background = vec3(0.);
