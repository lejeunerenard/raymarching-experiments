// cnoise loaded in frag.glsl & final-pass.glsl respectively

vec3 getBackground (in vec2 uv) {
  float coord = 1.0 * uv.y;

  return mix(vec3(0), #010101, coord);
}
vec3 background = vec3(0.);
