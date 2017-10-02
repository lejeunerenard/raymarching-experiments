#define PI 3.1415926536
#define TWO_PI 6.2831853072

precision highp float;

uniform vec2 resolution;
uniform float time;
uniform sampler2D base;
uniform sampler2D buffer;
uniform float wet;

#pragma glslify: import(./background)
#pragma glslify: cnoise2 = require(glsl-noise/classic/2d)

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

  vec2 UV = 2.0 * (uv - 0.5);
  float edgeWidth = 0.001;
  const float lineWidth = 0.002;
  float xCord = -UV.x - 0.85;
  float line = 0.0;

  for (int i = 0; i < 10; i++) {
    float localLine = smoothstep(-edgeWidth, 0., xCord);
    localLine *= smoothstep(edgeWidth + lineWidth, lineWidth, xCord);

    // End caps
    localLine *= smoothstep(0.8505, 0.85, abs(UV.y));

    float showMask = smoothstep(0.25, 0.251, cnoise2(vec2(451.34289 * xCord, time)));
    showMask += smoothstep(0.1, 0.11, cnoise2(vec2(0, 0.2 * time)));
    line += localLine * clamp(showMask, 0.0, 1.0);
    xCord += 0.025;
  }

  gl_FragColor = mix(gl_FragColor, vec4(#ffffff, 1.0), line);

  gl_FragColor.rgb = pow(gl_FragColor.rgb, gammaEnc);
}
