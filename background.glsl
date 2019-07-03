#pragma glslify: pModPolarBack = require(./hg_sdf/p-mod-polar-c.glsl)
#pragma glslify: rotMat2Back = require(./rotation-matrix2)

vec3 getBackground (in vec2 uv) {
  // Convert from [-1,1] -> [0, 1]
  vec2 coord = 0.5 * (uv.xy + vec2(1.0));

  vec3 color = vec3(0);

  const float r = 0.2 * 0.75;
  float n = 0.;

  // center
  vec2 q = uv;
  pModPolarBack(q, 6.);
  n = max(n, smoothstep(edge, 0., q.x - r));

  // b3
  q = uv - vec2(2. * r, 0);
  pModPolarBack(q, 6.);
  n = max(n, smoothstep(edge, 0., q.x - r));

  // b6
  q = uv - vec2(-2. * r, 0);
  pModPolarBack(q, 6.);
  n = max(n, smoothstep(edge, 0., q.x - r));

  // b1
  q = uv * rotMat2Back(0.33333 * PI);
  q -= vec2(2. * r, 0);
  pModPolarBack(q, 6.);
  n = max(n, smoothstep(edge, 0., q.x - r));

  // b4
  q = uv * rotMat2Back(0.33333 * PI);
  q -= vec2(-2. * r, 0);
  pModPolarBack(q, 6.);
  n = max(n, smoothstep(edge, 0., q.x - r));

  // b2
  q = uv * rotMat2Back(0.66667 * PI);
  q -= vec2(2. * r, 0);
  pModPolarBack(q, 6.);
  n = max(n, smoothstep(edge, 0., q.x - r));

  // b5
  q = uv * rotMat2Back(0.66667 * PI);
  q -= vec2(-2. * r, 0);
  pModPolarBack(q, 6.);
  n = max(n, smoothstep(edge, 0., q.x - r));

  color = vec3(n);

  color *= step(0.5, norT);

  return color;
}
vec3 background = vec3(0.);
