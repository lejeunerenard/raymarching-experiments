vec4 marchRef (in vec3 rayOrigin, in vec3 rayDirection, in float deltaT) {
  float t = 0.;
  float maxI = 0.;

  float trap = maxDistance;

  for (int i = 0; i < maxSteps / 3; i++) {
    vec3 d = map(rayOrigin + rayDirection * t, deltaT);
    if (d.x < epsilon) return vec4(t + d.x, d.y, float(i), d.z);
    t += d.x;
    maxI = float(i);
    trap = d.z;
    if (t > maxDistance) break;
  }
  return vec4(-1., 0., maxI, trap);
}
vec3 reflection (in vec3 ro, in vec3 rd, in float deltaT) {
  rd = normalize(rd);
  vec4 t = marchRef(ro + rd * .01, rd, deltaT);
  vec3 pos = ro + rd * t.x;
  vec3 color = vec3(0.);
  if (t.x > 0.) {
    vec3 nor = getNormal(pos, .0001, deltaT);
    color = diffuseColor(pos, nor, rd, t.y, t.w, deltaT);
    color /= max(1.0, pow(t.x, 1.0));
  } else {
    color = vec3(0.125); // 0.6 * vec3(getBackground(pos.xy));
  }

  return color;
}
vec3 reflection (in vec3 ro, in vec3 rd) {
  return reflection(ro, rd, 0.);
}

#pragma glslify: export(reflection)
