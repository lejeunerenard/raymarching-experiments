vec3 chen (in vec3 p, in float dt) {
  const float a = 40.0;
  const float b = 28.0;
  const float c = 3.0;

  float dxdt = dot(p.yx, vec2(a, -a));
  float dydt = (c - a) * p.x - p.x * p.z + c * p.y;
  float dzdt = p.x * p.y - b * p.z;

  return p + dt * vec3(dxdt, dydt, dzdt);
}

#pragma glslify: export(chen)
