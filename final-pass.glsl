#define PI 3.1415926536
#define TWO_PI 6.2831853072
#define saturate(x) clamp(x, 0.0, 1.0)

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
  color = mix(#00b3c8, color, 1. - l * .03125);
  // Dark
  color = mix(#9F6747, color, clamp(exp(l) * .755, 0., 1.));
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
  gl_FragColor.rgb = mix(gl_FragColor.rgb, colorBefore, 0.7);

  // Center origin
  uv = 2. * uv - 1.0;

  // Crop
  vec2 uvRainbow = abs(uv + vec2(0.5, 0)) * vec2(1, 0.45);
  float maskRainbow = 1.0 - smoothstep(0.4, 0.401, max(uvRainbow.x, uvRainbow.y));

  vec2 uvGrey = abs(uv - vec2(0.5, 0)) * vec2(1, 0.45);
  float maskGrey = 1.0 - smoothstep(0.4, 0.401, max(uvGrey.x, uvGrey.y));

  // White
  gl_FragColor.rgb += 10.0 * max(0., 1.0 - maskGrey - maskRainbow);

  gl_FragColor.rgb = saturate(gl_FragColor.rgb);

  // Gamma encode
  gl_FragColor.rgb = pow(gl_FragColor.rgb, gammaEnc);
}
