// source: https://www.shadertoy.com/view/ml2GWc

#define SHADOWS_QUALITY 1.0

// shadow function
float shadow(vec3 ro, vec3 rd, in float tmin, float tmax, in float time) {
  float t = tmin;
  for (float i=0.; i<100.; i++) {
    vec3 p = ro + rd*t;
    float h = map(p, time).x * SHADOWS_QUALITY;
    if (h<.001) return 0.;
    t += h;
    if (t >= tmax) break;
  }
  return 1.;
}

#pragma glslify: export(shadow)
