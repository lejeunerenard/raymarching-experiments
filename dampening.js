export function dampen (prev, next, factor) {
  return (1 - factor) * prev + factor * next
}
