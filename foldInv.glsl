void foldInv (inout vec2 p) {
  if (p.x - p.y < 0.) {
    float x1 = p.y;
    p.y = p.x;
    p.x = x1;
  }
}

#pragma glslify: export(foldInv)
