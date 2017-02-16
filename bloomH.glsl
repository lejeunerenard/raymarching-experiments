precision highp float;

uniform vec2 resolution;
uniform sampler2D buffer;

#pragma glslify: blur = require('glsl-fast-gaussian-blur')

void main() {
  vec2 uv = vec2(gl_FragCoord.xy / resolution.xy);
  gl_FragColor = blur(buffer, uv, resolution.xy, vec2(1., 0.));
}
