#define PI 3.1415926536
#define TWO_PI 6.2831853072
#define saturate(x) clamp(x, 0.0, 1.0)

precision highp float;

varying vec2 fragCoord;
uniform vec2 resolution;
uniform float time;
uniform sampler2D base;
uniform sampler2D buffer;
uniform sampler2D prevBuffer;
uniform float wet;

#pragma glslify: cnoise2 = require(glsl-noise/classic/2d)
#pragma glslify: import(./background)

void colorMap (inout vec3 color) {
  float l = length(vec4(color, 1.));
  // Light
  color = mix(#00b3c8, color, 1. - l * .05125);
  // Dark
  color = mix(#9F6747, color, clamp(exp(l) * .255, 0., 1.));
}

const float totalT = 8.0;
float modT = mod(time, totalT);
float norT = modT / totalT;
float cosT = TWO_PI / totalT * modT;

void main() {
  const vec3 gamma = vec3(2.2);
  const vec3 gammaEnc = vec3(0.454545);

  vec2 uv = vec2(gl_FragCoord.xy / resolution.xy);
  vec2 uvBackground = fragCoord.xy;
  background = getBackground(uvBackground);

  vec4 baseColor = texture2D(base, uv);
  baseColor.rgb = pow(baseColor.rgb, gamma);

  vec4 bufferColor  = texture2D(buffer, uv);
  bufferColor.rgb = pow(bufferColor.rgb, gamma);

  vec2 coord = uv; // 2.0 * (uv - 0.5);
  coord.y += 0.001;

  // coord *= 0.995;
  // coord = coord * 0.5 + 0.5;

  vec4 prevBufferColor = texture2D(prevBuffer, coord);
  prevBufferColor.rgb = pow(prevBufferColor.rgb, gamma);

  // Fade
  // prevBufferColor.a *= 0.9975;
  // prevBufferColor.a *= smoothstep(0., 0.001, prevBufferColor.a);

  vec4 result = mix(vec4(background, 1.), baseColor, baseColor.a);
  // vec4 result = mix(vec4(prevBufferColor.rgb, 1), baseColor, baseColor.a);

  gl_FragColor = wet * bufferColor + result;

  // gl_FragColor = mix(vec4(background, 1.), result, max(bufferColor.a, baseColor.a));
  // gl_FragColor = vec4(background + result.rgb, 1.);

  // Post process
  vec3 colorBefore = gl_FragColor.rgb;
  colorMap(gl_FragColor.rgb);
  gl_FragColor.rgb = mix(gl_FragColor.rgb, colorBefore, 0.85);

  // gl_FragColor.gb = pow(gl_FragColor.gb, vec2(1.5));

  // Gamma encode
  gl_FragColor.rgb = pow(gl_FragColor.rgb, gammaEnc);

  // 'Film' Noise
  gl_FragColor.rgb += 0.05 * (cnoise2(560. * uv + sin(uv + time) + 000.0 * vec2(time, 0.0)) + cnoise2(800. * uv + 253.5 * vec2(0., time)));
}
