precision highp float;

uniform vec2 resolution;
uniform sampler2D base;
uniform sampler2D buffer;

void main() {
  const float gamma = 1.7;

  vec2 uv = vec2(gl_FragCoord.xy / resolution.xy);
  vec4 result = texture2D(base, uv) + texture2D(buffer, uv);
  result.rgb = pow(result.rgb, vec3(1.0 / gamma));
  gl_FragColor = result;
}
