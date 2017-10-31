// cnoise loaded in frag.glsl & final-pass.glsl respectively

vec3 getBackground (in vec2 uv) {
  float coord = 1.0 * uv.y;
  return pow(#20080a, vec3(2.2));
}
vec3 background = vec3(0.);
