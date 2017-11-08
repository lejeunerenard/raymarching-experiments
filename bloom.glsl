precision highp float;

uniform vec2 resolution;
uniform vec2 direction;
uniform sampler2D buffer;
#define VARIABLE_BLOOM 1

#pragma glslify: blur = require(glsl-fast-gaussian-blur/13)
#ifdef VARIABLE_BLOOM
uniform float time;
#pragma glslify: blur9 = require(glsl-fast-gaussian-blur/9)
#pragma glslify: blur5 = require(glsl-fast-gaussian-blur/5)
#pragma glslify: cnoise2 = require(glsl-noise/classic/2d)
#pragma glslify: snoise2 = require(glsl-noise/simplex/2d)
#endif

void main() {
  vec2 uv = vec2(gl_FragCoord.xy / resolution.xy);
  #ifdef VARIABLE_BLOOM
  vec4 blurColor = vec4(0);

  float n = cnoise2(323.0 * uv + 10.0 * vec2(
    cnoise2(uv + vec2(time, 1.3 * time)),
    cnoise2(2.0 * uv + 2024.30 + vec2(time, 1.3 * time))));
  blurColor = blur(buffer, uv, resolution.xy, direction);
  blurColor = mix(blurColor, blur9(buffer, uv, resolution.xy, direction), smoothstep(0.5, 1.0, n));
  blurColor = mix(blurColor, blur5(buffer, uv, resolution.xy, direction), smoothstep(0.0, 0.5, n));

  gl_FragColor = blurColor;
  #else
  gl_FragColor = blur(buffer, uv, resolution.xy, direction);
  #endif
}
