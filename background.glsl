// cnoise loaded in frag.glsl & final-pass.glsl respectively

vec3 getBackground (in vec2 uv) {
  float coord = 1.0 * uv.y;

  return mix(#313030, #fff0f0, coord);
}
vec3 background = vec3(0.);
