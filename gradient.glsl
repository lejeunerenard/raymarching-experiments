vec3 gradientColors (in float v) {
  #define SIZE 5

  vec3 colors[SIZE];
  colors[0] = #E46BFF;
  colors[1] = #391B40;
  colors[2] = #AB51BF;
  colors[3] = #CD61E5;
  colors[4] = #72367F;

  v = mod(v, 1.);

  const float period = 1. / float(SIZE);

  vec3 color = colors[0];
  for (int i = 1; i < SIZE; i++) {
    float a = smoothstep(period * float(i - 1), period * float(i), v);
    color = mix(color, colors[i], a);
  }

  // Final loop around
  float a = smoothstep(period * float(SIZE - 1), 1., v);
  color = mix(color, colors[0], a);

  return color;
}
#pragma glslify: export(gradientColors)
