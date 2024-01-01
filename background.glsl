#pragma glslify: pModPolarBack = require(./hg_sdf/p-mod-polar-c.glsl)
#pragma glslify: pMod2 = require(./hg_sdf/p-mod2.glsl)

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
  // vec3 color = mix(vec3(0.0), vec3(0.05, 0.025, 0.05), bgIndex);
  // vec3 color = mix(vec3(0.1, 0.09, 0.09), 0.85 * vec3(0.1, 0.05, 0.15), bgIndex);
  // vec3 color = mix(vec3(0.7, 0.7, 0.9), vec3(1), bgIndex);
  // vec3 color = mix(0.9 * vec3(0, 1, 1), vec3(0, 0.4, 1), bgIndex);

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
  // vec3 color = mix(vec3(0.8, 0.9165, 1), vec3(0.5), 0.00);
  // color = mix(color, 0.7 * vec3(0, 1,0.8), bgIndex);

  // Solid colors
  // vec3 color = vec3(0.5);
  // vec3 color = vec3(0.7);
  // vec3 color = vec3(1, 0, 0);
  vec3 color = vec3(0);
  // vec3 color = mix(#5927F8, vec3(1), 0.20);
  // vec3 color = #FAABB9;

  // -- Patterns --

  // vec2 nQ = uv;
  // float warpScale = 2.;
  // nQ += warpScale * 0.10000 * cos( 2. * nQ.yx );
  // nQ += warpScale * 0.05000 * cos( 4. * nQ.yx );
  // nQ += warpScale * 0.02500 * cos( 8. * nQ.yx );

  // float n = dot(cos(23.78 * nQ), sin(71.29 * nQ));
  // n = dot(cos(vec2(0.2, 23.78) * nQ + n), sin(vec2(11.29, 0.7) * nQ));
  // n = dot(uv, vec2(0.1)) + 0.2 * n + 8. * dot(uv, uv * cos(vec2(1,3) * uv));

  // vec3 color = 0.5 + 0.5 * cos(TWO_PI * (n + vec3(uv, 1) + vec3(0, 0.33, 0.67)));

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
  // float period = 10.;
  // float n = sin(period * TWO_PI * dI);
  // // float cutoff = 0.8 * smoothstep(-0.5, 0.5, uv.y);
  // float cutoff = 0.0;
  // n = 1. - smoothstep(cutoff, cutoff + edge, n);
  // vec3 color = mix(vec3(0.0), vec3(1), n);

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
  // color *= 1.35; // Brighten

  // color = pow(color, vec3(2.2));

  // // Gradient tint
  // vec3 gradientColor = 0.5 + 0.5 * cos(TWO_PI * (vec3(uv, cos(dot(uv, vec2(1)))) + vec3(0, 0.33, 0.67)));
  // color *= mix(gradientColor, vec3(1), (0.70 + 0.2 * length(uv)));

  // color *= mix(vec3(0.95), vec3(1, 0.8, 0.8), saturate(1. - 0.95 * (uv.y + 1.0)));
  // color *= 1.2;
  // color = mix(color, vec3(0.5), 0.10); // desaturate
  // color = mix(color, vec3(1), mix(0.25, 0.80, 1. - bgIndex)); // lighten
  // color = mix(color, vec3(0), 0.3); // Darken
  // color *= 0.91;
  // color *= mix(vec3(1, 1.2, 1.2), vec3(1.3, 1.0, 1.0), bgIndex); // Tint
  color *= vec3(0.9, 0.8, 1.); // Tint

  return color;
}

vec3 getBackground (in vec2 uv) {
  return getBackground(uv, 0.);
}

vec3 background = vec3(0.);
