vec3 getBackground (in vec2 uv) {
  float coord = 1.0 * uv.y;

  return pow(mix(#FF998F, #FFB282, coord), vec3(2.2));
}
vec3 background = vec3(0.);
