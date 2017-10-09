// cnoise loaded in frag.glsl & final-pass.glsl respectively

vec3 getBackground (in vec2 uv) {
  float coord = 1.0 * uv.y;

  float n = 0.0;
  n += 0.2500 * cnoise2( 3.0 * uv);
  n += 0.1250 * cnoise2( 6.0 * uv);
  n += 0.0625 * cnoise2(12.0 * uv);

  coord += n;

  return mix(#000000, #030303, coord);
  return mix(#777777, #AAAAAA, coord);
}
vec3 background = vec3(0.);
