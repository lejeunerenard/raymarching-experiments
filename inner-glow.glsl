vec3 innerGlow (in float trap){
  float fGlow = clamp(trap * 0.1, 0.0, 1.0);
  fGlow = pow(fGlow, 3.5);

  return glowColor * 17.5 * fGlow;
}
