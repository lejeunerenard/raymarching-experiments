#pragma glslify: pModPolar = require(./hg_sdf/p-mod-polar-c.glsl)
#pragma glslify: rotMat2 = require(./rotation-matrix2)
const float backgroundR = 0.45;
float backgroundMask (in vec2 uv, in float r) {
  uv.y *= 2.0;
  uv *= rotMat2(PI * 0.25);
  vec2 absUV = abs(uv);
  return smoothstep(edge, 0., max(absUV.x, absUV.y) - r);
}

vec3 getBackground (in vec2 uv) {
  // Convert from [-1,1] -> [0, 1]
  vec2 coord = 0.5 * (uv.xy + vec2(1.0));

  float frameMask = backgroundMask(uv, backgroundR);
  vec2 backgroundUv = uv;
  backgroundUv.y *= 2.0;
  backgroundUv *= rotMat2(PI * 0.25);
  vec2 absUV = abs(backgroundUv);
  vec3 colorOffset = vec3(0, 0.33, 0.73);
  vec3 color = vec3(smoothstep(0.5 * edge, 0., abs(max(absUV.x, absUV.y) - backgroundR)));
  return color;
}
vec3 background = vec3(0.);
