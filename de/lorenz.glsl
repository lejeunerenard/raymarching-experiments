vec3 lorenz (in vec3 p, in float dt) {
  const float r = 28.0;
  const float s = 10.0;
  const float b = 2.666667; // 8/3

  float dxdt = s * dot(p.yx, vec2(1.0, -1.0));
  float dydt = p.x * ( r - p.z ) - p.y;
  float dzdt = p.x * p.y - b * p.z;

  return p + dt * vec3(dxdt, dydt, dzdt);
}

#pragma glslify: export(lorenz)
