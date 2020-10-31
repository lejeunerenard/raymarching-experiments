#pragma glslify: pModPolarBack = require(./hg_sdf/p-mod-polar-c.glsl)
#pragma glslify: pMod2 = require(./hg_sdf/p-mod2.glsl)
#pragma glslify: rotMat2Back = require(./rotation-matrix2)
// #pragma glslify: snoise2 = require(glsl-noise/simplex/2d)

vec3 getBackground (in vec2 uv, in float universe) {
  // Convert from [-1,1] -> [0, 1]
  vec2 coord = 0.5 * (uv.xy + vec2(1.0));

  // --- Gradient index ---
  // float bgIndex = saturate(0.5 * (uv.y + 1.0));
  float bgIndex = saturate(length(uv));
  // bgIndex = pow(bgIndex, 4.);
  // bgIndex = 1. - bgIndex; // Flip

  // --- Set colors / gradient ---
  // Gradients
  vec3 color = mix(vec3(0.05, 0.10, 0.2), vec3(0.0), bgIndex);
  // vec3 color = mix(vec3(0.05), vec3(0.0), bgIndex);
  // vec3 color = mix(vec3(1), 0.7 * #F8D261, bgIndex);

  // const vec3 bgColor = #F2900A;
  // vec3 color = mix(0.8 * bgColor, bgColor, bgIndex);
  // color = mix(color, vec3(1), 0.6);

  // Rainbow dark background color
  // float rainbowI = 0.5 * (atan(uv.y, uv.x) / PI + 1.);
  // rainbowI += 0.125 * sin(TWO_PI * rainbowI + sin(1. * TWO_PI * rainbowI));
  // rainbowI += norT;
  // float rainbowValue = 0.25 + 0.075 * cos(2. * TWO_PI * norT);
  // vec3 rainbow = 0.5 + 0.5 * cos(TWO_PI * (rainbowI + vec3(0, 0.33, 0.67)));
  // vec3 color = mix(rainbowValue * rainbow, vec3(0.0), bgIndex);

  // Hex color
  // vec3 color = mix(#FAC699, #FAC699, saturate(1.0 * bgIndex));

  // Solid colors
  // vec3 color = vec3(0);

  // -- Patterns --
  // // Grid
  // const float size = 0.015;
  // float sizeI = smoothstep(-1.0, 0.5, uv.y);
  // uv *= rotMat2Back(0.25 * PI);
  // vec2 c = pMod2(uv, vec2(size));
  // vec2 absQ = abs(uv);
  // float n = max(absQ.x, absQ.y) - (0.35 + 0.115 * sizeI) * size;
  // n = 1. - step(0.0, n);
  // vec3 color = 0.8 * vec3(n);

  // // Stripes
  // float dI = dot(uv, vec2(1));
  // float n = sin(47. * TWO_PI * dI);
  // float cutoff = 0.8 * smoothstep(-0.5, 0.5, uv.y);
  // n = 1. - smoothstep(cutoff, cutoff + edge, n);
  // vec3 color = 0.8 * vec3(n);

  // Manipulations
  // color = mix(color, #FFC070, saturate(smoothstep(0.0, 0.5, uv.y)));
  // color *= vec3(0.7, 0.75, 2.0); // Tint

  return color;
}

vec3 getBackground (in vec2 uv) {
  return getBackground(uv, 0.);
}

vec3 background = vec3(0.);
