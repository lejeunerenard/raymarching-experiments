#define PI 3.1415926536
#define TWO_PI 6.2831853072

precision highp float;

uniform vec2 resolution;
uniform float time;
uniform sampler2D base;
uniform sampler2D buffer;
uniform float wet;

#pragma glslify: import(./background)

void main() {
  const vec3 gamma = vec3(2.2);
  const vec3 gammaEnc = vec3(0.454545);

  vec2 uv = vec2(gl_FragCoord.xy / resolution.xy);
  background = getBackground(uv);

  vec4 baseColor = texture2D(base, uv);
  baseColor.rgb = pow(baseColor.rgb, gamma);

  vec4 bufferColor  = texture2D(buffer, uv);
  bufferColor.rgb = pow(bufferColor.rgb, gamma);

  vec4 result = wet * bufferColor + baseColor;

  // Gamma encode
  gl_FragColor = mix(vec4(background, 1.), result, baseColor.a);
  gl_FragColor.rgb = pow(gl_FragColor.rgb, gammaEnc);

  const float startT = 2.0;
  float mask = 1.0;

// Shot 1
  mask *= smoothstep(startT, startT + 1.0, time);
  mask *= smoothstep(startT + 7.0, startT + 5.0, time) + step(startT + 8.0, time);

// Shot 2
  mask *= (1.0 - step(startT + 7.99, time)) + smoothstep(startT + 8.1, startT + 9.0, time);
  mask *= smoothstep(startT + 15.0, startT + 14.0, time);

  gl_FragColor.rgb *= mask;
}
