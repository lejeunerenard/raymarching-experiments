#pragma glslify: pModPolarBack = require(./hg_sdf/p-mod-polar-c.glsl)
#pragma glslify: pMod2 = require(./hg_sdf/p-mod2.glsl)
#pragma glslify: rotMat2Back = require(./rotation-matrix2)
// #pragma glslify: snoise2 = require(glsl-noise/simplex/2d)

vec3 getBackground (in vec2 uv, in float universe) {
  // Convert from [-1,1] -> [0, 1]
  vec2 coord = 0.5 * (uv.xy + vec2(1.0));

  // --- Gradient index ---
  float bgIndex = saturate(0.5 * (uv.y + 1.0));
  // float bgIndex = saturate(length(uv));
  // bgIndex = pow(bgIndex, 4.);
  bgIndex = 1. - bgIndex; // Flip

  // --- Set colors / gradient ---
  // Gradients
  // vec3 color = mix(0.5 * vec3(0.3, 0.25, 0.1), vec3(0.0), bgIndex);
  // vec3 color = mix(vec3(0.25, 0.1, 0.1), vec3(0.15, 0.05, 0), bgIndex);
  vec3 color = mix(vec3(1.), vec3(0.9, 0.9, 0.975), bgIndex);
  // // Rainbow dark background color
  // float rainbowI = 0.5 * (atan(uv.y, uv.x) / PI + 1.);
  // rainbowI += 0.125 * sin(TWO_PI * rainbowI + sin(1. * TWO_PI * rainbowI));
  // rainbowI += norT;
  // float rainbowValue = 0.25 + 0.075 * cos(2. * TWO_PI * norT);
  // vec3 rainbow = 0.5 + 0.5 * cos(TWO_PI * (rainbowI + vec3(0, 0.33, 0.67)));
  // vec3 color = mix(rainbowValue * rainbow, vec3(0.0), bgIndex);
  // vec3 color = mix(#FF6A0A, #A11E0A, saturate(1.0 * bgIndex));
  // vec3 color = mix(#FF6A0A, #A11E0A, saturate(1.0 * bgIndex));

  // Solid colors
  // vec3 color = vec3(0);

  // Manipulations
  // color = mix(color, #FFC070, saturate(smoothstep(0.0, 0.5, uv.y)));
  // color *= vec3(0.7, 0.75, 2.0); // Tint

  return color;
}

vec3 getBackground (in vec2 uv) {
  return getBackground(uv, 0.);
}

vec3 background = vec3(0.);
