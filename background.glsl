#pragma glslify: pModPolarBack = require(./hg_sdf/p-mod-polar-c.glsl)
#pragma glslify: rotMat2Back = require(./rotation-matrix2)

vec3 getBackground (in vec2 uv) {
  // Convert from [-1,1] -> [0, 1]
  vec2 coord = 0.5 * (uv.xy + vec2(1.0));

  // vec3 color = mix(vec3(1.0), vec3(0.8), saturate(length(uv)));
  // vec3 color = mix(#FFC3A1, #FFD894, saturate(length(uv)));
  // vec3 color = mix(vec3(0.2), vec3(0.25), saturate(length(uv)));
  // vec3 color = mix(0.8 * vec3(0.9, 0, 1), vec3(0.0), mixI);
  // float n = dot(uv, vec2(37));
  // n = smoothstep(0., edge, sin(TWO_PI * n));
  // vec3 color = vec3(n);

  vec3 color = vec3(0.85); // vec3(0.25);

  return color;
}
vec3 background = vec3(0.);
