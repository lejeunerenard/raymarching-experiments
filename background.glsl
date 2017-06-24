vec3 getBackground (in vec2 uv) {
  float coord = 1.0 * uv.y;

  return pow(mix(#F7B3FF, #FFA6B7, coord), vec3(2.2));
}
vec3 background = vec3(0.);
