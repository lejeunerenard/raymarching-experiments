#pragma glslify: import(./generalized-polyhedra-include)

// p as usual, e exponent (p in the paper), r radius or something like that
float dodecahedral(vec3 p, float e, float r) {
  float s = pow(abs(dot(p,n14)),e);
  s += pow(abs(dot(p,n15)),e);
  s += pow(abs(dot(p,n16)),e);
  s += pow(abs(dot(p,n17)),e);
  s += pow(abs(dot(p,n18)),e);
  s += pow(abs(dot(p,n19)),e);
  s = pow(s, 1./e);
  return s-r;
}

#pragma glslify: export(dodecahedral)
