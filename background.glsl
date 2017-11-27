vec3 getBackground (in vec2 uv) {
  vec2 coord = 2.0 * (uv.xy - vec2(0.5));

  vec2 q = coord;
  float n = cnoise2(2.0 * q);
  q.x *= 2.5;
  n += 0.50 * cnoise2(4.0 * q);
  q.x *= 7.5;
  n += 0.25 * cnoise2(8.0 * q);

  return mix(vec3(0.001), #FF00FF, 0.007812 * n + 0.4 * saturate(1.0 - 0.98 * pow(length(coord), 0.0625)));
}
vec3 background = vec3(0.);
