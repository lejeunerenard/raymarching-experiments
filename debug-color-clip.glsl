void debugColor (inout vec3 color) {
  #ifdef enable
  if (color == vec3(0.)) {
    color = vec3(1., 0., 1.);
  }
  if (color == vec3(1.)) {
    color = vec3(0., 1., 0.);
  }
  #endif
}

#pragma glslify: export(debugColor)
