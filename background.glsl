#pragma glslify: pModPolarBack = require(./hg_sdf/p-mod-polar-c.glsl)
#pragma glslify: pMod2 = require(./hg_sdf/p-mod2.glsl)
#pragma glslify: rotMat2Back = require(./rotation-matrix2)
// #pragma glslify: snoise2 = require(glsl-noise/simplex/2d)

vec3 getBackground (in vec2 uv, in float universe) {
  // Convert from [-1,1] -> [0, 1]
  vec2 coord = 0.5 * (uv.xy + vec2(1.0));

  // --- Gradient index ---
  // float bgIndex = saturate(0.5 * (1. * uv.y + 1.0));
  float bgIndex = saturate(length(uv));
  // bgIndex = pow(bgIndex, 4.);
  // bgIndex += 0.1 * dot(sin(6. * uv), vec2(1));
  bgIndex = 1. - bgIndex; // Flip

  // Metadata

  // --- Set colors / gradient ---
  // Gradients
  // vec3 color = mix(#501E5B, vec3(0.0), bgIndex);
  vec3 color = mix(vec3(0.10), vec3(0.05, 0, 0), bgIndex);
  // vec3 color = mix(vec3(0.9), vec3(0.6), bgIndex);
  // vec3 color = 1.2 * mix(vec3(0.3), vec3(0.90), bgIndex);
  // vec3 color = mix(#BD96E5, #9498F2, bgIndex);
  // color *= mix(colors1, vec3(1), 1. - length(coord));
  // color = mix(color, vec3(1), 0.30);

  // const vec3 bgColor = #F2900A;
  // vec3 color = mix(0.8 * bgColor, bgColor, bgIndex);
  // color = mix(color, vec3(1), 0.6);

  // // Rainbow dark background color
  // float rainbowI = 0.5 * (atan(uv.y, uv.x) / PI + 1.);
  // rainbowI += 0.125 * sin(TWO_PI * rainbowI + sin(1. * TWO_PI * rainbowI));
  // rainbowI += norT;
  // float rainbowValue = 1.;
  // vec3 rainbow = 0.5 + 0.5 * cos(TWO_PI * (rainbowI + vec3(0, 0.33, 0.67)));
  // vec3 color = mix(rainbowValue * rainbow, vec3(0.0), bgIndex);

  // Hex color
  // vec3 color = mix(#34A0E0, #60E083, saturate(1. * bgIndex));
  // vec3 color = mix(0.25 * #91acff, 0.5 * #a19cff, saturate(1. * bgIndex));
  // color *= 1.5;
  // vec3 color = 0.25 * #91acff;

  // Solid colors
  // vec3 color = vec3(0.5);
  // vec3 color = vec3(0, 0, 1);
  // vec3 color = vec3(0);
  // vec3 color = mix(#F82759, vec3(0.5), 0.30);

  // -- Patterns --

  // // Grid
  // const float size = 0.007;
  // float sizeI = smoothstep(-1.0, 0.25, uv.y);
  // uv *= rotMat2Back(0.35 * PI);
  // vec2 c = pMod2(uv, vec2(size));
  // vec2 absQ = abs(uv);
  // float n = max(absQ.x, absQ.y) - (0.35 + 0.115 * sizeI) * size;
  // n = 1. - step(0.0, n);
  // vec3 color = 0.8 * vec3(n);

  // // Stripes
  // vec2 axis = vec2(1);
  // float dI = dot(uv, axis);
  // float period = 32.;
  // float n = sin(period * TWO_PI * dI);
  // // float cutoff = 0.8 * smoothstep(-0.5, 0.5, uv.y);
  // float cutoff = 0.;
  // n = 1. - smoothstep(cutoff, cutoff + edge, n);
  // vec3 color = mix(vec3(0.0), vec3(1), n); // vec3(1.00 * n); // mix(#FAC011, #001FAD, n);

  // // Dots
  // float size = 0.0265;
  // vec2 c = pMod2(uv, vec2(size));
  // float dotN = length(uv) - 0.3 * size;
  // // n = sin(54. * TWO_PI * dI);
  // // float dotCutoff = 0.8 * smoothstep(-0.5, 0.5, uv.y);
  // float dotCutoff = 0.;
  // dotN = smoothstep(dotCutoff, dotCutoff + edge, dotN);
  // vec3 color = vec3(1.00 * dotN);
  // // color = mix(mix(primeColor, vec3(1), dotN), color, m);

  // Manipulations
  // color = mix(color, #FFC070, saturate(smoothstep(0.0, 0.5, uv.y)));
  // color *= vec3(0.7, 4.0, 1.3); // Tint
  // color *= 1.35; // Brighten

  // color = pow(color, vec3(2.2));

  color *= vec3(0.79, 1.0, 0.8);
  // color = mix(color, vec3(1), 0.175);

  return color;
}

vec3 getBackground (in vec2 uv) {
  return getBackground(uv, 0.);
}

vec3 background = vec3(0.);
