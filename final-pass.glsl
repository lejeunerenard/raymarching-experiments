#define PI 3.1415926536
#define TWO_PI 6.2831853072

precision highp float;

uniform vec2 resolution;
uniform float time;
uniform sampler2D base;
uniform sampler2D buffer;
uniform float wet;

#pragma glslify: cnoise2 = require(glsl-noise/classic/2d)
#pragma glslify: import(./background)

void colorMap (inout vec3 color) {
  float l = length(vec4(color, 1.));
  // Light
  color = mix(#043210, color, 1. - l * .0625);
  // Dark
  color = mix(#ef78FF, color, clamp(exp(l) * .325, 0., 1.));
}

void main() {
  const vec3 gamma = vec3(2.2);
  const vec3 gammaEnc = vec3(0.454545);

  vec2 uv = vec2(gl_FragCoord.xy / resolution.xy);
  background = getBackground(uv);

  vec4 baseColor = texture2D(base, uv);
  baseColor.rgb = pow(baseColor.rgb, gamma);

  vec4 bufferColor  = texture2D(buffer, uv);
  bufferColor.rgb = pow(bufferColor.rgb, gamma);

  vec4 result = mix(vec4(background, 1.), baseColor, baseColor.a);
  gl_FragColor = wet * bufferColor + result;

  // gl_FragColor = mix(vec4(background, 1.), result, max(bufferColor.a, baseColor.a));
  // gl_FragColor = vec4(background + result.rgb, 1.);

  // Post process
  vec3 colorBefore = gl_FragColor.rgb;
  colorMap(gl_FragColor.rgb);
  gl_FragColor.rgb = mix(gl_FragColor.rgb, colorBefore, 0.5);

  // Gamma encode
  gl_FragColor.rgb = pow(gl_FragColor.rgb, gammaEnc);
}
