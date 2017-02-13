void foldNd (inout vec3 z, vec3 n1) {
  z-=2.0 * min(0.0, dot(z, n1)) * n1;
}

#pragma glslify: export(foldNd)
