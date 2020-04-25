#pragma glslify: pModPolarBack = require(./hg_sdf/p-mod-polar-c.glsl)
#pragma glslify: rotMat2Back = require(./rotation-matrix2)

vec3 getBackground (in vec2 uv) {
  // Convert from [-1,1] -> [0, 1]
  vec2 coord = 0.5 * (uv.xy + vec2(1.0));

  // vec3 color = mix(vec3(0.7), vec3(1.0), saturate(length(uv)));


  const vec3 darkColor = 0.5 * #6C71CC;
  vec3 color = mix(#FB94FF, darkColor, 1. - saturate(0.4 * uv.y));
  vec3 innerColor = mix(#FB94FF, 0.6 * #6C71CC, 1. - saturate(-0.7 * uv.y));

  vec2 absUV = vec2(1, 0.45) * abs(uv);
  float maxD = max(absUV.x, absUV.y) - 0.225;
  float innerColorI = smoothstep(edge, 0., maxD);
  color = mix(color, innerColor, saturate(innerColorI));

  // vec3 color = mix(vec3(0.25), vec3(0.0), saturate(length(uv)));
  // vec3 color = mix(0.8 * vec3(0.9, 0, 1), vec3(0.0), mixI);

  // vec3 color = 0.5 * #6C71CC;
  // vec3 color = vec3(0.125);

  return color;
}
vec3 background = vec3(0.);
