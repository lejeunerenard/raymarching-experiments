float hue2IORPolynomial (in float hue, in float greenIOR, in float n1, in float null) {
  float x = 0.0005556 * hue - 0.25002; // 1.0 / (360.0 * 2.0) * ( 0.4 * hue - 180.0);
  float between = 0.05 * (x + 4.0) * (x + 2.0) * (x + 1.0); // * (x - 1.0); // * (x - 3.0) + 2.0;
  float ior = n1/(greenIOR + between);
  return ior;
}

#pragma glslify: export(hue2IORPolynomial)
