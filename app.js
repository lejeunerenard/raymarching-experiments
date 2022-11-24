import createShader from 'gl-shader'
import createTexture from 'gl-texture2d'
import createFBO from 'gl-fbo'

import ndarray from 'ndarray'
import TWEEN from 'tween.js'
import makeContext from 'gl-context'
import { rot4 } from './utils'
import drawTriangle from 'a-big-triangle'

import defined from 'defined'
import assert from 'assert'
import { vec3, mat4 } from 'gl-matrix'
import { getLuminance, getColorWFixedLuminance } from './luminance'

const dpr = 1.0 / window.devicePixelRatio

const TWO_PI = 2 * Math.PI
// const PHI = (1 + Math.sqrt(5)) / 2

const MANDELBOX = false
const BLOOM = true
const BLOOM_PASSES = 2
const BLOOM_WET = 1.0
const BLOOM_MIN_BRIGHTNESS = 0.99

// Initialize shell
export default class App {
  constructor (options = {}) {
    let canvas = document.createElement('canvas')
    document.body.appendChild(canvas)
    canvas.style.display = 'none'

    let gl = makeContext(canvas, { preserveDrawingBuffer: true })

    // enable extensions
    var ext = gl.getExtension('OES_standard_derivatives')
    if (!ext) {
      throw new Error('derivatives not supported')
    }

    this.LOOKAT = true

    this.presets = {}
    const preset = {
      offset: {
        x: 1.513,
        y: -0.702,
        z: 1.292
      },
      d: 0.52,
      scale: 1.3087,
      rot2angle: [1.18, 0, 5.902],
      cameraAngles: [0.102, 0.032, 0.025]
    }

    this.d = preset.d
    this.cameraRo = vec3.fromValues(0.22, 0.8, 1.76)
    this.offsetC = [0.339, -0.592, 0.228, 0.008]

    this.colors1 = [168, 141, 198]
    this.colors2 = [75, 24, 17]
    // this.getEqualLuminance(this.colors1, this.colors2, 0)

    // Ray Marching Parameters
    this.epsilon = preset.epsilon || 0.000001

    // Fractal parameters
    this.offset = (preset.offset)
      ? vec3.fromValues(preset.offset.x, preset.offset.y, preset.offset.z)
      : vec3.fromValues(0, 0, 0)
    this.scale = preset.scale
    this.rot2angle = preset.rot2angle || [0, 0, 0]
    this.cameraAngles = preset.cameraAngles || [0, 0, 0]

    this.angle1C = -0.117
    this.angle2C = -1.4645
    this.angle3C = 2.103

    // this.setupAnimation(preset)

    this.glInit(gl)

    // Audio
    this.audioFFT = 32
    this.audioTexArray = new Uint8Array(1 * this.audioFFT)
    this.audioNday = ndarray(this.audioTexArray, [this.audioFFT, 1])
    this.audioTex = createTexture(gl, this.audioNday)
    this.pulseGoal = 0

    // Capturing state
    this.capturing = defined(options.capturing, false)

    this.loaded = Promise.resolve()
      .then(() => {
        this.setupAudio()
      })

    // Scene Rendering
    this.sceneRenderer = options.sceneRenderer

    Object.assign(this, {
      canvas,
      gl
    })
  }

  getEqualLuminance (reference, second, position) {
    second[position] = null
    second[position] = getColorWFixedLuminance(
      second[0], second[1], second[2],
      getLuminance(reference))
  }

  getDimensions () {
    let width = this.width || window.innerWidth
    let height = this.height || window.innerHeight
    return [dpr * width, dpr * height]
  }

  setupFBOs (gl) {
    let dim = this.getDimensions()
    this.state = [
      createFBO(gl, dim, { depth: false }),
      createFBO(gl, dim, { depth: false }),
      createFBO(gl, dim, { depth: false }),
      createFBO(gl, dim, { depth: false }),
      createFBO(gl, dim, { depth: false }) ]

    this.state[0].color.magFilter = gl.LINEAR
    this.state[0].color.minFilter = gl.LINEAR
    this.state[1].color.magFilter = gl.LINEAR
    this.state[1].color.minFilter = gl.LINEAR
    this.state[2].color.magFilter = gl.LINEAR
    this.state[2].color.minFilter = gl.LINEAR
    this.state[3].color.magFilter = gl.LINEAR
    this.state[3].color.minFilter = gl.LINEAR
    this.state[4].color.magFilter = gl.LINEAR
    this.state[4].color.minFilter = gl.LINEAR
  }

  setupAnimation (preset) {
    let self = this
    // Epsilon Animation
    let eps1 = new TWEEN.Tween(this)
    eps1
      .to({ epsilon: 0.001 }, 10 * 1000)
      .easing(TWEEN.Easing.Quadratic.InOut)
    // eps1.start(0)

    // Camera location animation
    let ob = {
      x: self.cameraRo[0],
      y: self.cameraRo[1],
      z: self.cameraRo[2]
    }
    function updatePos () {
      self.cameraRo[0] = this.x
      self.cameraRo[1] = this.y
      self.cameraRo[2] = this.z
    }

    const totalTime = 10

    let cameraPosTween = new TWEEN.Tween(ob)
    cameraPosTween
      .delay(2.5 * 1000)
      .to({ x: 0, y: 0, z: self.cameraRo[2] }, 2.5 * 1000)
      .easing(TWEEN.Easing.Quadratic.Out)
      .onUpdate(updatePos)

    let cameraPosTween2 = new TWEEN.Tween(ob)
    cameraPosTween2
      .delay(2.5 * 1000)
      .to({ x: self.cameraRo[0], y: self.cameraRo[1], z: self.cameraRo[2] }, 2.5 * 1000)
      .easing(TWEEN.Easing.Quadratic.Out)
      .onUpdate(updatePos)

    cameraPosTween.chain(cameraPosTween2)
    cameraPosTween2.chain(cameraPosTween)
    // cameraPosTween.start(0)

    // Camera rotation
    function updateRot () {
      self.cameraAngles[0] = this[0]
      self.cameraAngles[1] = this[1]
      self.cameraAngles[2] = this[2]
    }

    let rotObj = [...this.cameraAngles]
    let camRotTween1 = new TWEEN.Tween(rotObj)
    camRotTween1
      .to([-0.328, this.cameraAngles[1], this.cameraAngles[2]], 4 * 1000)
      .onUpdate(updateRot)
      .easing(TWEEN.Easing.Linear.None)
    let camRotTween2 = new TWEEN.Tween(rotObj)
    camRotTween2
      .to([...this.cameraAngles], 4 * 1000)
      .onUpdate(updateRot)
      .easing(TWEEN.Easing.Linear.None)

    camRotTween1.chain(camRotTween2)
    camRotTween2.chain(camRotTween1)
    // camRotTween1.start(0)

    // Animation Fractal
    let rotTween1 = new TWEEN.Tween(this.rot2angle)
    rotTween1
      .delay(1 * 1000)
      .to([this.rot2angle[0], 2.255, this.rot2angle[2]], 7 * 1000)
      .easing(TWEEN.Easing.Quadratic.InOut)
    let rotTween2 = new TWEEN.Tween(this.rot2angle)
    rotTween2
      .delay(1 * 1000)
      .to([...this.rot2angle], 6 * 1000)
      .easing(TWEEN.Easing.Quadratic.InOut)

    rotTween1.chain(rotTween2)
    rotTween2.chain(rotTween1)

    rotTween1.start(0)

    // Scale Tween
    let scaleTween1 = new TWEEN.Tween(this)
    scaleTween1
      .delay(0 * 1000)
      .to({ scale: 1.4911 }, 5 * 1000)
      .easing(TWEEN.Easing.Quadratic.InOut)
    let scaleTween2 = new TWEEN.Tween(this)
    scaleTween2
      .delay(0 * 1000)
      .to({ scale: this.scale }, 5 * 1000)
      .easing(TWEEN.Easing.Quadratic.InOut)

    scaleTween1.chain(scaleTween2)
    scaleTween2.chain(scaleTween1)
    // scaleTween1.start(0)

    // Offset Tween
    let offsetTween1 = new TWEEN.Tween(this.offset)
    offsetTween1
      .to([ this.offset[0], this.offset[1], Math.PI ], 30 * 1000)
      .easing(TWEEN.Easing.Quadratic.InOut)
    let offsetTween2 = new TWEEN.Tween(this.offset)
    offsetTween2
      .to([ 1.993, this.offset[1], -0.654 ], 5 * 1000)
      .easing(TWEEN.Easing.Quadratic.InOut)
    let offsetTween3 = new TWEEN.Tween(this.offset)
    offsetTween3
      .to([ ...this.offset ], 5 * 1000)
      .easing(TWEEN.Easing.Quadratic.InOut)

    // offsetTween1.chain(offsetTween2)
    // offsetTween2.chain(offsetTween3)
    // offsetTween3.chain(offsetTween1)

    // offsetTween1.start(0)

    // Angle1C Tween
    let angle1CTween1 = new TWEEN.Tween(this)
    angle1CTween1
      .to({ angle1C: 0.5166 }, 20 * 1000)
      .easing(TWEEN.Easing.Linear.None)
    // angle1CTween1.start(0)

    // // Angle2C Tween
    // let angle2CTween1 = new TWEEN.Tween(this)
    // angle2CTween1
    //   .to({ angle2C: TWO_PI * 2 }, 15 * 1000)
    //   .easing(TWEEN.Easing.Quadratic.InOut)

    // angle2CTween1.start(0)

    // Angle3C Tween
    // let angle3CTween1 = new TWEEN.Tween(this)
    // angle3CTween1
    //   .to({ angle3C: 1.25 }, 20 * 1000)
    //   .easing(TWEEN.Easing.Quadratic.InOut)

    // angle3CTween1.start(0)
  }

  setupAudio () {
    return new Promise((resolve, reject) => {
      const audioCtx = new (window.AudioContext || window.webkitAudioContext)()

      const output = audioCtx.createGain()
      output.gain.setValueAtTime(0.2, audioCtx.currentTime)
      output.connect(audioCtx.destination)

      this.analyser = audioCtx.createAnalyser()
      this.analyser.fftSize = this.audioFFT
    })
  }

  enableEvents () {
    if (this.capturing) return
    this.resizeBound = this.resizeBound || this.resize.bind(this)
    window.addEventListener('resize', this.resizeBound, true)
    window.addEventListener('vrdisplaypresentchange', this.resizeBound, true)
  }

  disposeEvents () {
    window.removeEventListener('resize', this.resizeBound, true)
    window.removeEventListener('vrdisplaypresentchange', this.resizeBound, true)
  }

  setupShader (property, shader, gl) {
    let showError = (e) => {
      let pre = document.createElement('pre')
      pre.classList.add('glsl-error')

      let code = document.createElement('code')
      code.innerHTML = e.message
      pre.appendChild(code)
      document.body.appendChild(pre)
    }

    try {
      this[property] = createShader(gl, shader.vertex, shader.fragment)
      shader.on('change', () => {
        // Remove existing errors
        const errors = document.body.querySelectorAll('.glsl-error')
        for (const error of errors) {
          error.parentNode.removeChild(error)
        }

        try {
          this[property] = createShader(gl, shader.vertex, shader.fragment)
        } catch (e) {
          if (e.name === 'GLError') {
            showError(e)
          } else {
            throw e
          }
        }
      })
    } catch (e) {
      if (e.name === 'GLError') {
        showError(e)
      } else {
        throw e
      }
    }
  }

  glInit (gl) {
    // Turn off depth test
    gl.disable(gl.DEPTH_TEST)

    // Create fragment shader
    this.setupShader('shader', require('./shaders/frag.shader'), gl)

    this.setupShader('bright', require('./shaders/bright.shader'), gl)
    this.setupShader('bloom', require('./shaders/bloom.shader'), gl)
    this.setupShader('finalPass', require('./shaders/final-pass.shader'), gl)

    this.currentState = 1
    this.setupFBOs(gl)

    this.shader.attributes.position.location = 0
  }

  kifsM (t = 0, scale = this.scale, offset = this.offset) {
    this.shader.uniforms.scale = scale
    this.shader.uniforms.offset = offset

    // Scale and Offset
    let _kifsM

    if (MANDELBOX) {
      _kifsM = mat4.fromValues(
        1, 0, 0, -offset[0],
        0, 1, 0, -offset[1],
        0, 0, 1, -offset[2],
        0, 0, 0, 1)
    } else {
      _kifsM = mat4.fromValues(
        scale, 0, 0, -offset[0] * (scale - 1),
        0, scale, 0, -offset[1] * (scale - 1),
        0, 0, scale, -offset[2] * (scale - 1),
        0, 0, 0, 1)
    }

    const angleX = this.rot2angle[0]
    this.shader.uniforms.rot = angleX
    const axisX = vec3.fromValues(1, 0, 0)
    mat4.multiply(_kifsM, rot4(axisX, angleX), _kifsM)

    // Y-centric
    const angleY = this.rot2angle[1]
    const axisY = vec3.fromValues(0, 1, 0)
    mat4.multiply(_kifsM, rot4(axisY, angleY), _kifsM)

    // Z-centric
    const angleZ = this.rot2angle[2]
    const axisZ = vec3.fromValues(0, 0, 1)
    mat4.multiply(_kifsM, rot4(axisZ, angleZ), _kifsM)

    return _kifsM
  }

  resize (e) {
    let canvas = this.canvas
    let dim = this.getDimensions()
    canvas.width = dim[0]
    canvas.height = dim[1]
    canvas.style.width = dim[0] + 'px'
    canvas.style.height = dim[1] + 'px'

    this.state[0].shape = dim
    this.state[1].shape = dim
    this.state[2].shape = dim
    this.state[3].shape = dim
    this.state[4].shape = dim
  }

  tick (t) {
    this.shader.bind()

    this.update(t)
    this.render(t)
  }

  getCamera (t) {
    t /= 1000
    let cameraMatrix = mat4.create()

    // LookAt
    if (this.LOOKAT) {
      mat4.lookAt(cameraMatrix, this.cameraRo, vec3.fromValues(0, 0, 0), vec3.fromValues(0, 1, 0))
    } else {
      const angleX = this.cameraAngles[0]
      const axisX = vec3.fromValues(1, 0, 0)
      mat4.multiply(cameraMatrix, rot4(axisX, angleX), cameraMatrix)

      // Y-centric
      const angleY = this.cameraAngles[1]
      const axisY = vec3.fromValues(0, 1, 0)
      mat4.multiply(cameraMatrix, rot4(axisY, angleY), cameraMatrix)

      // Z-centric
      const angleZ = this.cameraAngles[2]
      const axisZ = vec3.fromValues(0, 0, 1)
      mat4.multiply(cameraMatrix, rot4(axisZ, angleZ), cameraMatrix)
    }

    this.cameraMatrix = cameraMatrix
    return [this.cameraRo, cameraMatrix]
  }

  update (t) {
    t = (window.time !== undefined) ? window.time : t
    TWEEN.update(t)

    this.shader.uniforms.epsilon = this.epsilon

    let updates = this.getCamera(t)
    this.shader.uniforms.cameraRo = updates[0]
    this.shader.uniforms.cameraMatrix = (updates[1])

    this.shader.uniforms.kifsM = this.kifsM(t, this.scale, this.offset)
    this.shader.uniforms.offsetC = this.offsetC

    this.shader.uniforms.angle1C = this.angle1C
    this.shader.uniforms.angle2C = this.angle2C
    this.shader.uniforms.angle3C = this.angle3C

    this.shader.uniforms.colors1 = [this.colors1[0] / 255, this.colors1[1] / 255, this.colors1[2] / 255]
    // Update colors2 based on colors1 luminance
    // this.getEqualLuminance(this.colors1, this.colors2, 0)
    this.shader.uniforms.colors2 = [this.colors2[0] / 255, this.colors2[1] / 255, this.colors2[2] / 255]

    this.shader.uniforms.d = this.d
  }

  bloomBlur (gl, t) {
    let dim = this.getDimensions()

    // Brightness pass
    let base = this.state[0].color[0]
    this.state[1].bind()
    this.bright.bind()
    this.bright.uniforms.minBright = BLOOM_MIN_BRIGHTNESS
    this.bright.uniforms.buffer = base.bind(0)
    this.bright.uniforms.resolution = dim
    drawTriangle(gl)

    for (let i = 0; i < BLOOM_PASSES; i++) {
      // Horizontal Blur
      let brightLayer = this.state[1].color[0]
      this.state[2].bind()

      this.bloom.bind()
      this.bloom.uniforms.buffer = brightLayer.bind(1)
      this.bloom.uniforms.resolution = dim
      this.bloom.uniforms.direction = [1, 0]
      this.bloom.uniforms.time = this.getTime(t)
      drawTriangle(gl)

      // Vertical Blur
      let prev = this.state[2].color[0]
      this.state[1].bind()

      this.bloom.uniforms.buffer = prev.bind(2)
      this.bloom.uniforms.resolution = dim
      this.bloom.uniforms.direction = [0, 1]
      this.bloom.uniforms.time = this.getTime(t)
      drawTriangle(gl)
    }

    // Additive blending
    this.currentState = (this.currentState) ? 0 : 1
    this.state[3 + this.currentState].bind()
    this.finalPass.bind()
    this.finalPass.uniforms.base = base.bind(3)
    this.finalPass.uniforms.buffer = this.state[1].color[0].bind(4)
    this.finalPass.uniforms.prevBuffer = this.state[3 + ((this.currentState + 1) % 2)].color[0].bind(5)
    this.finalPass.uniforms.resolution = dim
    this.finalPass.uniforms.time = this.getTime(t)
    this.finalPass.uniforms.wet = BLOOM_WET
    this.finalPass.uniforms.colors1 = [this.colors1[0] / 255, this.colors1[1] / 255, this.colors1[2] / 255]
    this.finalPass.uniforms.colors2 = [this.colors2[0] / 255, this.colors2[1] / 255, this.colors2[2] / 255]
    drawTriangle(gl)

    // Render again as framebuffer
    gl.bindFramebuffer(gl.FRAMEBUFFER, null)
    drawTriangle(gl)
  }

  getTime (t) {
    return window.time || t / 1000
  }

  render (t) {
    let { shader, gl } = this

    if (BLOOM) {
      this.state[0].bind()
    }

    shader.uniforms.time = this.getTime(t)
    shader.uniforms.BLOOM = BLOOM
    this.sceneRenderer(shader, t)

    if (BLOOM) {
      this.bloomBlur(gl, t)
    }
  }

  run () {
    assert(this.sceneRenderer, 'A sceneRenderer is required')

    this.canvas.style.display = null
    this.resize()

    this.enableEvents()
  }

  stop () {
    this.canvas.style.display = 'none'
    this.disposeEvents()
  }
}
