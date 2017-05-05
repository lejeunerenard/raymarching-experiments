vec3 quartz (in vec3 p, in float angle) {
  vec3 quartzBody = vec3(sdHexPrism(p.xzy - vec3(0., -0.30, 0.0), vec2(0.5, 3.)), 1., 0.);

  const float dodecaScale = 1.00;
  vec3 dp = p;
  dp.y *= 0.416;
  dp *= rotationMatrix(vec3(0.000,PHI,1.), PI / 5. + angle) / dodecaScale;
  vec3 d = vec3(dodeca(dp, 1.) * dodecaScale, 1., 0.);

  return dMax(quartzBody, d);
  // return d;
}

#pragma glslify: export(quartz)
