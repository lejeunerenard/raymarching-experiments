#pragma glslify: pModPolarBack = require(./hg_sdf/p-mod-polar-c.glsl)
#pragma glslify: rotMat2Back = require(./rotation-matrix2)

vec3 getBackground (in vec2 uv) {
  // Convert from [-1,1] -> [0, 1]
  vec2 coord = 0.5 * (uv.xy + vec2(1.0));

  // vec3 color = mix(vec3(1), vec3(0.8), saturate(length(uv)));
  // vec3 color = mix(vec3(0.5), vec3(0.7), saturate(length(uv)));
  vec3 color = mix(vec3(0.), vec3(0.2), saturate(length(uv)));

  float stripeSpeed = 4.;
  float n = sin(TWO_PI * stripeSpeed * 8. * dot(fragCoord.xy, vec2(0, 1)));
  n = smoothstep(edge, 0., n);

  float invN = 0.;
  // float invN = sin(TWO_PI * stripeSpeed * 8. * dot(fragCoord.xy, vec2(1, 0)));
  // invN = smoothstep(edge, 0., invN);

  n = mix(n, 1. - n, invN);

  color = vec3(n);

  return color;
}
vec3 background = vec3(0.);
