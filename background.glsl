// cnoise loaded in frag.glsl & final-pass.glsl respectively

vec3 getBackground (in vec2 uv) {
  float coord = 1.0 * uv.y;

  return vec3(0.01);
}
vec3 background = vec3(0.);
