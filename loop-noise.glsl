float loopNoise (in vec3 q, in float startNoise, in float endNoise) {
  const float timeScale = 1.;
  // q must be normalized on `.x`
  float n1 = noise(vec3(timeScale, 1, 1) * q);
  float n2 = noise(vec3(timeScale, 0, 0) + vec3(-timeScale, 1, 1) * q);
  float noiseIndex = clamp((q.x - startNoise) / (endNoise - startNoise), 0., 1.);
  return mix(n1, n2, noiseIndex);
}

#pragma glslify: export(loopNoise)
