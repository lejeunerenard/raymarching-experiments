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

  // Transitions
  const float startTime = 0.0;

  // Scene 1
  gl_FragColor.rgb *= smoothstep(startTime + 0.0,  startTime + 2.0, time); // Fade In
  gl_FragColor.rgb *= 1.0 - smoothstep(startTime + 4.0, startTime + 5.0, time) // Fade Out
    + step(startTime + 5.1, time); // Disable at 5.1s

  // Scene 2
  gl_FragColor.rgb *= smoothstep(startTime + 5.2, startTime + 8.0, time) // Fade In
   + 1.0 - step(startTime + 5.1, time); // Enable at 5.1s
  gl_FragColor.rgb += step(startTime + 14.0, startTime + time) * 2.0 * (0.5 + 0.5 * sin(TWO_PI * 0.2 * pow(time - 5.0, 3.0)));
  gl_FragColor.rgb *= 1.0 - step(startTime + 16.5, startTime + time);
}
