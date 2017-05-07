precision highp float;

uniform vec2 resolution;
uniform sampler2D base;
uniform sampler2D buffer;

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

  vec4 result = baseColor + bufferColor;

  // Gamma encode
  gl_FragColor = mix(vec4(background, 1.), result, baseColor.a);
  gl_FragColor.rgb = pow(result.rgb, gammaEnc);
}
