// Needs:
// - map(vec2 pos, vec2 size)
// - float numberOfNeighbors
// - float maxDistance

#pragma glslify: pMod2 = require(../hg_sdf/p-mod2.glsl)

float neighborGrid (in vec2 q, in vec2 size) {
  vec2 c = pMod2(q, size);

  float n = maxDistance;

  for (float x = -numberOfNeighbors; x < numberOfNeighbors + 1.; x++)
  for (float y = -numberOfNeighbors; y < numberOfNeighbors + 1.; y++) {
    vec2 shift = vec2(x, y);
    float cell = map(q - size * shift, c + shift);
    n = min(n, cell);
  }

  return n;
}

#pragma glslify: export(neighborGrid)
