float loopNoise (in vec3 q, in float startNoise, in float endNoise) {
  // q must be normalized on `.x`
  float n1 = noise(q);
  float n2 = noise(vec3(1, 0, 0) + vec3(-1, 1, 1) * q);
  float noiseIndex = clamp((q.x - startNoise) / (endNoise - startNoise), 0., 1.);
  return mix(n1, n2, noiseIndex);
}

#pragma glslify: export(loopNoise)
