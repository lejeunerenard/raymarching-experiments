import assert from 'assert'
import drawTriangle from 'a-big-triangle'

class DefaultSceneRenderer {
  constructor (gl) {
    assert(gl, 'A gl context is required')
    this.gl = gl
    this.width = window.innerWidth
    this.height = window.innerHeight
  }

  resize (width, height) {
    this.width = width
    this.height = height
  }

  render (shader, t) {
    let { width, height, gl } = this
    shader.uniforms.resolution = [width, height]
    gl.viewport(0, 0, width, height)

    var projectionMatrix = this.projectionMatrix(80, 1, 1)
    shader.uniforms.projectionMatrix = projectionMatrix

    drawTriangle(gl)
  }

  projectionMatrix (fov, aspect, zoom) {
    const near = 0.01
    const far = 10000
    const DEG2RAD = Math.PI / 180

    const skew = 0
    const filmGauge = 35

    var top = near * Math.tan(
    DEG2RAD * 0.5 * fov) / zoom

    var height = 2 * top
    var width = aspect * height
    var left = -0.5 * width

    const filmWidth = filmGauge * Math.min(aspect, 1)
    if (skew !== 0) left += near * skew / filmWidth

    return this.makeFrustum(left, left + width, top - height, top, near, far)
  }

  makeFrustum (left, right, bottom, top, near, far) {
    var te = []
    var x = 2 * near / (right - left)
    var y = 2 * near / (top - bottom)

    var a = (right + left) / (right - left)
    var b = (top + bottom) / (top - bottom)
    var c = -(far + near) / (far - near)
    var d = -2 * far * near / (far - near)

    te[0] = x
    te[4] = 0
    te[8] = a
    te[12] = 0

    te[1] = 0
    te[5] = y
    te[9] = b
    te[13] = 0

    te[2] = 0
    te[6] = 0
    te[10] = c
    te[14] = d

    te[3] = 0
    te[7] = 0
    te[11] = -1
    te[15] = 0

    return te
  }
}

module.exports = DefaultSceneRenderer
