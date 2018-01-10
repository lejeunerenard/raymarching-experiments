float hue2IORExp (in float hue, in float greenIOR, in float n1, in float between) {
  float relPos = 0.00089 * hue - 0.5; // 0.5 / 360.0 * (hue - 180.0);
  float ior = n1/(greenIOR + 0.50 * sign(relPos) * pow(3.5, abs(relPos)));
  return ior;
}

#pragma glslify: export(hue2IORExp)
