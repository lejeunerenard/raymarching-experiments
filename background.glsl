vec3 getBackground (in vec2 uv) {
  float coord = 1.0 - 1.0 * uv.y;
  coord *= coord * coord;

  return hsv(vec3(coord * 0.125, 1.0, 1.0));
}
vec3 background = vec3(0.);
