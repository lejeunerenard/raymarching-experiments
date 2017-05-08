#pragma glslify: import(./generalized-polyhedra-include)

// p as usual, e exponent (p in the paper), r radius or something like that
float icosahedral(vec3 p, float e, float r) {
  float s = pow(abs(dot(p,n4)),e);
  s += pow(abs(dot(p,n5)),e);
  s += pow(abs(dot(p,n6)),e);
  s += pow(abs(dot(p,n7)),e);
  s += pow(abs(dot(p,n8)),e);
  s += pow(abs(dot(p,n9)),e);
  s += pow(abs(dot(p,n10)),e);
  s += pow(abs(dot(p,n11)),e);
  s += pow(abs(dot(p,n12)),e);
  s += pow(abs(dot(p,n13)),e);
  s = pow(s, 1./e);
  return s-r;
}

#pragma glslify: export(icosahedral)
