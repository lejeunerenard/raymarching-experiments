// cnoise loaded in frag.glsl & final-pass.glsl respectively

vec3 getBackground (in vec2 uv) {
  float coord = 1.0 * uv.y;

  float n = 0.0;
  n += 0.500 * cnoise2( 3.0 * uv);
  n += 0.250 * cnoise2( 6.0 * uv);
  n += 0.125 * cnoise2(12.0 * uv);

  coord += 0.75 * n;

  return mix(#E89694, #FFAFD0, coord);
}
vec3 background = vec3(0.);
