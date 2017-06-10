float hue2IORExp (in float hue, in float greenIOR, in float n1, in float between) {

  float relPos = 0.002778 * hue - 0.5; // 1.0 / 360.0 * (hue - 180.0);
  float ior = n1/(greenIOR + 0.3 * sign(relPos) * pow(3.0, abs(relPos)));
  return ior;
}

#pragma glslify: export(hue2IORExp)
