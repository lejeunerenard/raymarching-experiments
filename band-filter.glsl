float bandFilter (in float t, in float start, in float end) {
  return smoothstep(start, end, t) - step(end, t);
}
#pragma glslify: export(bandFilter)
