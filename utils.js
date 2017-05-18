import { mat4, vec3 } from 'gl-matrix'

export function isIOS () {
  return !!navigator.platform && /iPad|iPhone|iPod/.test(navigator.platform)
}
export function isAndroid () {
  let ua = navigator.userAgent.toLowerCase()
  return ua.indexOf('android') > -1
}

export function rot4 (axis, angle) {
  const s = Math.sin(angle)
  const c = Math.cos(angle)
  const oc = 1.0 - c
  vec3.normalize(axis, axis)

  return mat4.fromValues(
    oc * axis[0] * axis[0] + c, oc * axis[0] * axis[1] - axis[2] * s, oc * axis[2] * axis[0] + axis[1] * s, 0.0,
    oc * axis[0] * axis[1] + axis[2] * s, oc * axis[1] * axis[1] + c, oc * axis[1] * axis[2] - axis[0] * s, 0.0,
    oc * axis[2] * axis[0] - axis[1] * s, oc * axis[1] * axis[2] + axis[0] * s, oc * axis[2] * axis[2] + c, 0.0,
    0.0, 0.0, 0.0, 1.0)
}
