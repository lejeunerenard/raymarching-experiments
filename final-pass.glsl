precision highp float;

uniform vec2 resolution;
uniform sampler2D base;
uniform sampler2D buffer;

void main() {
  const float gamma = 1.5;

  vec2 uv = vec2(gl_FragCoord.xy / resolution.xy);
  vec4 baseColor = texture2D(base, uv);
  vec4 result = baseColor + texture2D(buffer, uv);
  result.rgb = pow(result.rgb, vec3(1.0 / gamma));
  gl_FragColor = mix(vec4(#eeeefe, 1.), result, baseColor.a);
}
