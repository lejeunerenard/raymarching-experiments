// cnoise loaded in frag.glsl & final-pass.glsl respectively

vec3 getBackground (in vec2 uv) {
  float coord = 1.0 * uv.y;

  float v = sin(TWO_PI * (dot(uv, vec2(9, -9)) + 0.2 * time));
  v = smoothstep(0.0, 0.15, v);
  return vec3(v);
}
vec3 background = vec3(0.);
