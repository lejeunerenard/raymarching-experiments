vec3 cluster (in vec3 x) {
  const float scale = 0.4;
  x /= scale;

  const float angleSpan = PI * 0.85;

  vec3 res = vec3(1000.0, 1., 0.);
  for (int k=-1; k<=1; k++) {
    for (int j=-1; j<=1; j++) {
      for (int i=-1; i<=1; i++) {
        vec3 b = 1.5 * vec3(i, j, k);

        vec3 p = x + b;
        p *= rotationMatrix(vec3(
          noise(b),
          noise(1.1 * b - 20.3),
          noise(0.9 * b + 414.9)), angleSpan * noise(20.0 * b));

        vec3 d = quartz(p, PI / 4.0 * noise(13.0 * b));
        d.x *= scale;
        res = dMin( res, d );
      }
    }
  }
  return res;
}
#pragma glslify: export(cluser)
