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

uniform vec3 colors1;
uniform vec3 colors2;

#pragma glslify: import(./time)
const float edge = 0.0024;

#pragma glslify: cnoise2 = require(glsl-noise/classic/2d)
#pragma glslify: import(./background)

void colorMap (inout vec3 color) {
  float l = length(vec4(color, 1.));
  // Light
  color = mix(#00b3c8, color, 1. - l * .05125);
  // Dark
  color = mix(#9F6747, color, clamp(exp(l) * .255, 0., 1.));
}

void main() {
  modT = mod(time, totalT);
  norT = modT / totalT;
  cosT = TWO_PI / totalT * modT;
  const vec3 gamma = vec3(2.2);
  const vec3 gammaEnc = vec3(0.454545);

  vec2 uv = vec2(gl_FragCoord.xy / resolution.xy);
  vec2 uvBackground = fragCoord.xy;
  background = getBackground(uvBackground, step(0.5, norT));

  vec4 baseColor = texture2D(base, uv);
  baseColor.rgb = pow(baseColor.rgb, gamma);

  vec4 bufferColor  = texture2D(buffer, uv);
  bufferColor.rgb = pow(bufferColor.rgb, gamma);

  // uv = uvBackground;
  vec2 coord = uv - 0.5; // 2.0 * (uv - 0.5);
  vec2 norm = coord;
  norm *= 0.950;
  coord = norm + 0.5;
  // coord -= 0.001 * sin(vec2(cosT) + vec2(0, PI * 0.5));

  vec4 prevBufferColor = texture2D(prevBuffer, coord);
  prevBufferColor.rgb = pow(prevBufferColor.rgb, gamma);

  // Fade
  // prevBufferColor.a *= 0.9975;
  // prevBufferColor.a *= smoothstep(0., 0.001, prevBufferColor.a);

  vec4 result = mix(vec4(background, 1.), baseColor, baseColor.a);
  // vec4 result = mix(prevBufferColor, baseColor, baseColor.a);

  gl_FragColor = wet * bufferColor + result;
  // gl_FragColor = vec4(vec3(baseColor.a), 1);
  // if (norT > 0.2) {
  //   gl_FragColor = vec4(prevBufferColor.rgb, 1);
  // }

  // gl_FragColor = vec4(coord, 0, 1);
  // if (abs(coord.x - 0.5) < 0.01) {
  //   gl_FragColor = vec4(1, 0, 1, 1);
  // }
  // if (abs(coord.y - 0.5) < 0.01) {
  //   gl_FragColor = vec4(0, 1, 1, 1);
  // }

  // gl_FragColor = vec4( (0.995 * (coord - 0.5)) + 0.5, 0, 1 );

  // gl_FragColor = mix(vec4(background, 1.), result, max(bufferColor.a, baseColor.a));
  // gl_FragColor = vec4(background + result.rgb, 1.);

  // Post process
  vec3 colorBefore = gl_FragColor.rgb;
  colorMap(gl_FragColor.rgb);
  gl_FragColor.rgb = mix(gl_FragColor.rgb, colorBefore, 0.85);

  // gl_FragColor.gb = pow(gl_FragColor.gb, vec2(1.1));

  // Gamma encode
  // gl_FragColor.rgb = pow(gl_FragColor.rgb, gammaEnc);

  // 'Film' Noise
  // const float filmNoiseScale = 1.5;
  // gl_FragColor.rgb += 0.030 * (cnoise2(filmNoiseScale * 560. * uv + sin(uv + time)) + cnoise2(filmNoiseScale * 800. * uv + 253.5 * vec2(0., time)));
}
