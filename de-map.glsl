vec3 de (in vec3 x, in float endt, in float dt) {
  float t = 0.0;

  for (int i = 0; i < MAX_IT; i ++) {
    if (t >= endt) break;

    x = field(x, dt);

    t += dt;
  }

  return x;
}

#pragma glslify: export(de)
