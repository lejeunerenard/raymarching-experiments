vec3 getBackground (in vec2 uv) {
  // Convert from [-1,1] -> [0, 1]
  vec2 coord = 0.5 * (uv.xy + vec2(1.0));

  vec3 color = mix(vec3(0.0), vec3(0.06) , coord.x);

  float circ1 = length(uv + vec2(0.05, -0.05)) - 0.65;
  vec3 circ1Color = mix(pow(#FF6BE6, vec3(2.2)), pow(#52FFF5, vec3(2.2)), coord.x);
  circ1Color += 0.5 * (0.5 + 0.5 * cos(TWO_PI * (norT + vec3(0, 0.33, 0.67))));
  color = mix(color, circ1Color, smoothstep(edge, 0., circ1));

  float circ2 = length(uv - vec2(0.3, -0.3)) - 0.35;
  vec3 circ2Color = mix(pow(#FFE882, vec3(2.2)), pow(#52FFF5, vec3(2.2)), coord.x - 0.3);
  circ2Color += 0.5 * (0.5 + 0.5 * cos(TWO_PI * (norT + vec3(0, 0.1, 0.2))));
  color = mix(color, circ2Color, smoothstep(edge, 0., circ2));

  color = pow(color, vec3(1.5));

  return color;
}
vec3 background = vec3(0.);
