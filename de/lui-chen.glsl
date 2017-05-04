vec3 equation (in vec3 p, in float dt) {
  const float a = 2.4;
  const float b = -3.78;
  const float c = 14.0;
  const float d = -11.0;
  const float e = 4.0;
  const float f = 5.58;
  const float r = 1.0;

  float dxdt = dot(p.yx, vec2(a, b)) + c * p.y * p.z;
  float dydt = dot(p.yz, vec2(d, -1)) + e * p.x * p.z;
  float dzdt = f * p.z + r * p.x * p.y;

  return p + dt * vec3(dxdt, dydt, dzdt);
}

#pragma glslify: export(equation)
