#pragma glslify: import(./generalized-polyhedra-include)

// p as usual, e exponent (p in the paper), r radius or something like that
float icosahedral(vec3 p, float e, float r) {
  // Notes: shows and is super pointy. screen capture is at ~33fps
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

  // dot product equivalent (or so i think)
  // note: didnt show anything big nor small that I could see
  // float s = abs(dot(p,n4));
  // s += abs(dot(p,n5));
  // s += abs(dot(p,n6));
  // s += abs(dot(p,n7));
  // s += abs(dot(p,n8));
  // s += abs(dot(p,n9));
  // s += abs(dot(p,n10));
  // s += abs(dot(p,n11));
  // s += abs(dot(p,n12));
  // s += abs(dot(p,n13));
  // s = pow(s, 1./e);

  // // Maxi metric version
  // // Notes: shows and is super pointy. screen capture is at ~35fps
  // float s = abs(dot(p,n4));
  // s = max(s, abs(dot(p,n5)));
  // s = max(s, abs(dot(p,n6)));
  // s = max(s, abs(dot(p,n7)));
  // s = max(s, abs(dot(p,n8)));
  // s = max(s, abs(dot(p,n9)));
  // s = max(s, abs(dot(p,n10)));
  // s = max(s, abs(dot(p,n11)));
  // s = max(s, abs(dot(p,n12)));
  // s = max(s, abs(dot(p,n13)));

  return s-r;
}

#pragma glslify: export(icosahedral)
