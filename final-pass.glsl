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

  // Post Graphics
  vec2 UV = 2.0 * (uv - 0.5);
  float line = smoothstep(0.026, 0.025, abs(dot(UV, vec2(-1, 1))));
  line *= smoothstep(1.505, 1.5, abs(dot(UV, vec2(1))));
  gl_FragColor = mix(gl_FragColor, vec4(#ffffff, 1.0), line);
  // gl_FragColor = vec4(vec3(max(UV.x, UV.y)), 1.0);

  gl_FragColor.rgb = pow(gl_FragColor.rgb, gammaEnc);
}
