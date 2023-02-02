float gyroidTriangle (in vec3 p, in float thickness) {
  float gyroid = dot(triangleWave(p), triangleWave(p.yzx + 0.5));
  return abs(gyroid) - thickness;
}

#pragma glslify: export(gyroidTriangle)
