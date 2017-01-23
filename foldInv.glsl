void foldInv (inout vec2 p) {
  if (p.x - p.y < 0.) {
    float x1 = p.y;
    p.y = p.x;
    p.x = x1;
  }
}
void bfold (inout vec2 p) {
  const vec2 n = vec2(1., -1.);
  p -= min(0., dot(p, n)) * n;
}

#pragma glslify: export(bfold)
