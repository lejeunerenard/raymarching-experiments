#ifndef Iterations
#define Iterations 15
#endif

const vec3 un = vec3(1., -1., 0.);

vec2 aexion( inout vec3 p ) {
  float trap = maxDistance;

  const float delta = 4.;
  vec4 CT = vec4(
    abs(dot(p, un.xxx) - delta),
    abs(dot(p, un.yyx) - delta),
    abs(dot(p, un.yxy) - delta),
    abs(dot(p, un.xyy) - delta)
  );
  vec4 V = vec4(0.);
  float V2 = 0.;
  float dr = 2.;

  for (int i = 0; i < Iterations; i++) {
    V = clamp(V, -1., 1.) * 2. - V;
    V2 = dot(V,V);

    float c = clamp(max(.25 / V2, .25), 0., 1.) * 4.;
    V *= c;
    dr /= c;

    V = V * 2. + CT;
    dr *= .5;

    trap = min(trap, length(p));

    if (V2 > 3600.) break;
  }

  return vec2(sqrt(V2) * dr, trap);
}

#pragma glslify: export(aexion)
