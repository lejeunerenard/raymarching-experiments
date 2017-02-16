precision highp float;

uniform vec2 resolution;
uniform sampler2D buffer;
uniform float minBright;

void main() {
  vec2 uv = vec2(gl_FragCoord.xy / resolution.xy);

  vec4 color = texture2D(buffer, uv);
  // Check whether fragment output is higher than threshold, if so output as brightness color
  float brightness = dot(color.rgb, vec3(0.2126, 0.7152, 0.0722));
  if (brightness > .5) {
    gl_FragColor = vec4(color.rgb, 1.0);
  } else {
    gl_FragColor = vec4(vec3(0.), 1.);
  }
}
