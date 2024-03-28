float gyroid (in vec3 p, in float thickness) {
  float gyroid = dot(sin(p), cos(p.yzx));
  return abs(gyroid) - thickness;
}

#pragma glslify: export(gyroid)
