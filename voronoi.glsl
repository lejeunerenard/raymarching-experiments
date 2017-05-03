vec2 hash( vec2 p  ) {
  p = vec2(dot(p,vec2(127.1,311.7)),dot(p,vec2(269.5,183.3)));
  return fract(sin(p)*18.5453);
}

vec3 hash( vec3 p  ) {
  p = vec3(
    dot(p,vec3(127.1,311.7, 131.)),
    dot(p,vec3(269.5,183.3, 1.394)),
    dot(p,vec3(137.1,728.1, 91.230))
  );
  return fract(sin(p)*18.5453);
}

float voronoi(in vec2 x) {
  vec2 p = floor(x);
  vec2 f = fract(x);

  float res = 8.0;
  for (int j=-1; j<=1; j++) {
    for (int i=-1; i<=1; i++) {
      vec2 b = vec2(i, j);
      vec2  r = vec2(b) - f + hash(p + b);
      float d = dot( r, r );

      res = min( res, d );
    }
  }
  return sqrt( res );
}

float voronoi(in vec3 x) {
  vec3 p = floor(x);
  vec3 f = fract(x);

  float res = 8.0;
  for (int k=-1; k<=1; k++) {
    for (int j=-1; j<=1; j++) {
      for (int i=-1; i<=1; i++) {
        vec3 b = vec3(i, j, k);
        vec3 r = abs(vec3(b) - f + hash(p + b));
        float d = max(max(r.x, r.y), r.z);

        res = min( res, d );
      }
    }
  }
  return res;
}

#pragma glslify: export(voronoi)
