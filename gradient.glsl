vec3 gradientColors (in float v) {
  #define SIZE 3

  vec3 colors[SIZE];
  colors[0] = #AD27FF;
  colors[1] = #2C66FF;
  colors[2] = #39FFE0;

  v = mod(v, 1.);

  const float period = 1. / float(SIZE);

  vec3 color = colors[0];
  for (int i = 1; i < SIZE; i++) {
    float a =  smoothstep(period * float(i - 1), period * float(i), v);
    color = mix(color, colors[i], a);
  }

  // Final loop around
  float a =  smoothstep(period * float(SIZE - 1), 1., v);
  color = mix(color, colors[0], a);

  return color;
}
#pragma glslify: export(gradientColors)
