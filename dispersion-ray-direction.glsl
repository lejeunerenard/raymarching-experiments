float chunkedHueIOR (in float hue, in float greenIOR, in float n1, in float between) {
  float redIORRatio = n1/(greenIOR - 1.0 * between);
  float yellowIORRatio = n1/(greenIOR - 0.5 * between);
  float greenIORRatio = n1/greenIOR;
  float cyanIORRatio = n1/(greenIOR + 0.5 * between);
  float blueIORRatio = n1/(greenIOR + 1.0 * between);
  float purpleIORRatio = n1/(greenIOR + 1.5 * between);

  float ior = mix(redIORRatio, yellowIORRatio, smoothstep(0.0, 60.0, hue));
  ior = mix(ior, greenIORRatio, smoothstep(60.0, 120.0, hue));
  ior = mix(ior, cyanIORRatio, smoothstep(120.0, 180.0, hue));
  ior = mix(ior, blueIORRatio, smoothstep(180.0, 240.0, hue));
  ior = mix(ior, purpleIORRatio, smoothstep(240.0, 270.0, hue));

  return ior;
}

#pragma glslify: export(chunkedHueIOR)
