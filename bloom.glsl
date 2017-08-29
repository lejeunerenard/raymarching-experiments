precision highp float;

uniform vec2 resolution;
uniform vec2 direction;
uniform sampler2D buffer;

#pragma glslify: blur = require(glsl-fast-gaussian-blur/13)

void main() {
  vec2 uv = vec2(gl_FragCoord.xy / resolution.xy);
  gl_FragColor = blur(buffer, uv, resolution.xy, direction);
}
