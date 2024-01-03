// Needs:
// - map(vec2 pos, vec2 id)
// - float numberOfNeighbors
// - float maxDistance

#pragma glslify: pMod2 = require(../hg_sdf/p-mod2.glsl)

vec2 neighborGrid (in vec2 q, in vec2 size) {
  vec2 c = pMod2(q, size);

  float n = maxDistance;
  float m = -1.; // default material

  for (float x = -numberOfNeighbors; x < numberOfNeighbors + 1.; x++)
  for (float y = -numberOfNeighbors; y < numberOfNeighbors + 1.; y++) {
    vec2 shift = vec2(x, y);
    vec2 cell = map(q - size * shift, c + shift);
    if (cell.x < n) {
      m = cell.y;
    }
    n = min(n, cell.x);
  }

  return vec2(n, m);
}

// vec2 neighborGrid (in vec3 q, in vec2 size) {
//   vec2 c = pMod2(q.xz, size);

//   float n = maxDistance;
//   float m = -1.; // default material

//   for (float x = -numberOfNeighbors; x < numberOfNeighbors + 1.; x++)
//   for (float y = -numberOfNeighbors; y < numberOfNeighbors + 1.; y++) {
//     vec2 shift = vec2(x, y);
//     vec2 spaceShift = size * shift;
//     vec2 cell = map(q - vec3(spaceShift.x, 0, spaceShift.y), c + shift);
//     if (cell.x < n) {
//       m = cell.y;
//     }
//     n = min(n, cell.x);
//   }

//   return vec2(n, m);
// }

#pragma glslify: export(neighborGrid)
