precision highp float;

uniform vec2 resolution;
uniform sampler2D base;
uniform sampler2D buffer;

#pragma glslify: blur = require(glsl-fast-gaussian-blur/13)

void main() {
  vec2 uv = vec2(gl_FragCoord.xy / resolution.xy);
  vec4 bloom = blur(buffer, uv, resolution.xy, vec2(0., 1.));

  gl_FragColor = texture2D(base, uv) + bloom;
}
