float sigmoid ( in float x ) {
  const float L = 1.0;
  const float k = 1.0;
  const float x0 = 4.0;

  x *= 8.0; // Scale so x [0, 1]

  return L / ( 1.0 + exp(-k * (x - x0)) );
}

float hue2IORSigmoid (in float hue, in float greenIOR, in float n1, in float between) {
  float relPos = sigmoid(hue * 0.002778);
  float ior = n1/(greenIOR + relPos * between);
  return ior;
}

#pragma glslify: export(hue2IORSigmoid)
