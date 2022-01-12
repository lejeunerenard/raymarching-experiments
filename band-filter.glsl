float bandFilter (in float t, in float start, in float end) {
  float mid = (end - start) * 0.5;
  mid = min(mid, 0.001);
  return smoothstep(start, start + mid, t) - smoothstep(end - mid, end, t);
}
#pragma glslify: export(bandFilter)
