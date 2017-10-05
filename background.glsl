// cnoise loaded in frag.glsl & final-pass.glsl respectively

vec3 getBackground (in vec2 uv) {
  float coord = 1.0 * uv.y;

  vec2 uvN = 7.0 * uv;
  float n = 0.0;
  n += 0.25 * cnoise2(uvN);
  uvN *= 2.01;
  n += 0.125 * cnoise2(uvN);
  uvN *= 2.11;
  n += 0.0625 * cnoise2(uvN);
  uvN *= 2.03;

  coord += 0.5 * n;

  return mix(#777777, #EEEEEE, coord);
}
vec3 background = vec3(0.);
