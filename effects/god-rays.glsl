#ifdef saturate
// saturate is being scoped to the module so alias
#define saturate_0(x) saturate(x)
#else
#define saturate(x) clamp(x, 0.0, 1.0)
#endif

vec4 godRays ( in vec3 rayOrigin, in vec3 rayDirection, in vec4 t, in vec2 uv, in vec4 color, in float generalT) {
  const float maxGodRaySteps = 2.;

  // TODO try with gPos
  vec3 pos = rayOrigin + rayDirection * t.x;

  float godRayGlowIntensity = 1.;
  float godRayMask = 1.; // saturate(pow(1. - dot(uv, uv), 3. * (1. - godRayGlowIntensity)));
  if (godRayMask == 0.) {
    return color;
  }

  // Replace this with t.y or something
  float distanceToPos = distance(rayOrigin, pos);
  vec3 beamStep = rayDirection * distanceToPos / maxGodRaySteps;

  float illumination = 0.;
  float jitter = 0.0; // 0.05 * noise(uv);
  for (float i = 0.; i < maxGodRaySteps; i++) {
    vec3 samplePoint = rayOrigin + beamStep * (i + jitter); // Jitter so smooth sampling point.
    float shadow = softshadow(samplePoint, normalize(sunPos), 0.01, 1.0, generalT);
    // shadow = max(0.1, shadow);
    // illumination += saturate(shadow); //  / (distance(rayOrigin, sunPos));
    illumination += shadow; //  / (distance(rayOrigin, sunPos));
  }

  // illumination /= sqrt(maxGodRaySteps);
  // illumination /= pow(maxGodRaySteps, 0.85);
  // illumination *= 0.2 + 0.8 * godRayGlowIntensity;
  illumination *= godRayMask;

  float existingAlpha = color.a;
  // illumination = pow(illumination, 1. - 0.3 * (1. - existingAlpha));
  // color = mix(color, vec4(sunColor, 1.), illumination);
  // color = vec4(vec3(illumination), 1.);
  // color += vec4(sunColor, 1.) * illumination;

  // Debug : Show only illumination value
  color = vec4(sunColor * illumination, 1.);

  return color;
}

#pragma glslify: export(godRays)
