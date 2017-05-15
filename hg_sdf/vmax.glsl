// Maximum/minumum elements of a vector
float vmax(vec2 v) {
  return max(v.x, v.y);
}

float vmax(vec3 v) {
  return max(max(v.x, v.y), v.z);
}

float vmax(vec4 v) {
  return max(max(v.x, v.y), max(v.z, v.w));
}

#pragma glslify: export(vmax)
