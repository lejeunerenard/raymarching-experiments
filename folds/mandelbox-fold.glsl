vec4 mandelboxFold (in vec4 z, inout float minD) {
  vec4 z0 = vec4(z.xyz, 1.);

  float minRadius2 = minRadius * minRadius;
  vec4 scalevec = vec4(s, s, s, abs(s)) / minRadius2;

  for (int i = 0; i < trap; i++) {
    // Box Fold
    z.xyz = clamp(z.xyz, -foldLimit, foldLimit) * 2. - z.xyz;

    // Ball fold
    float r2 = dot(z.xyz, z.xyz);
    z.xyzw *= clamp(max(minRadius2/r2, minRadius2), 0., 1.);

    z *= rotM;

    z.xyzw = z*scalevec + z0;

    minD = min(minD, trapCalc(z));
  }

  return z;
}

#pragma glslify: export(mandelboxFold)
