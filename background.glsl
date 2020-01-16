#pragma glslify: pModPolarBack = require(./hg_sdf/p-mod-polar-c.glsl)
#pragma glslify: rotMat2Back = require(./rotation-matrix2)

vec3 getBackground (in vec2 uv) {
  // Convert from [-1,1] -> [0, 1]
  vec2 coord = 0.5 * (uv.xy + vec2(1.0));

  vec3 color = vec3(0); // mix(vec3(1), vec3(0.6), saturate(length(uv)));

  const float size = 0.1;
  coord = uv;
  vec2 c = floor((coord + size*0.5)/size);
  coord = mod(coord + size*0.5,size) - size*0.5;

  vec2 absCoord = abs(coord);
  float n = abs(min(absCoord.x, absCoord.y)) - size * 0.00625;

  n = 1. - step(0., n);

  n *= 1. - step(0., max(absCoord.x, absCoord.y) - size * 0.2);

  color = vec3(n);

  return color;
}
vec3 background = vec3(0.);
