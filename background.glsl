// cnoise loaded in frag.glsl & final-pass.glsl respectively

vec3 getBackground (in vec2 uv) {
  float coord = 1.0 * uv.y;

  return mix(vec3(0.1, 0.101, 0.101), vec3(0.91, 0.9, 0.9), coord);
}
vec3 background = vec3(0.);
