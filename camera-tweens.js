import TWEEN from 'tween.js'

export function cameraOrbit (ro, radius, from, to, period) {
  let tweened = from.slice()
  let tween = new TWEEN.Tween(tweened)
  tween
    .to(to, 7000)
    .onUpdate(function () {
      ro[0] = radius * Math.sin(this[0])
      ro[1] = radius * Math.cos(this[1])
    })
    .easing(TWEEN.Easing.Quadratic.InOut)
    .onComplete(function (object) {
      this[0] = from[0]
      this[1] = from[1]
    })
  return tween
}
