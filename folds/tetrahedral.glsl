vec3 fold (in vec3 q) {
  if (dot(q.xy, vec2(1, -1)) < 0.) q.yx = vec2(1) * q.xy;
  if (dot(q.xz, vec2(1, -1)) < 0.) q.zx = vec2(1) * q.xz;
  if (dot(q.yz, vec2(1, -1)) < 0.) q.zy = vec2(1) * q.yz;

  if (dot(q.xy, vec2(1)) < 0.) q.yx = vec2(-1) * q.xy;
  if (dot(q.xz, vec2(1)) < 0.) q.zx = vec2(-1) * q.xz;
  if (dot(q.yz, vec2(1)) < 0.) q.zy = vec2(-1) * q.yz;

  return q;
}

#pragma glslify: export(fold)
