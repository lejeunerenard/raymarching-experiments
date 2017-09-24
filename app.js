const glslify = require('glslify')

import createShader from 'gl-shader'
import createTexture from 'gl-texture2d'
import createFBO from 'gl-fbo'

import ndarray from 'ndarray'
import fit from 'canvas-fit'
import TWEEN from 'tween.js'
import makeContext from 'gl-context'
import { rot4 } from './utils'
import SoundCloud from 'soundcloud-badge'
// import Analyser from 'web-audio-analyser'
import drawTriangle from 'a-big-triangle'

import defined from 'defined'
import { vec3, mat4 } from 'gl-matrix'

const dpr = 0.75 * Math.min(2, defined(window.devicePixelRatio, 1))
const CLIENT_ID = 'ded451c6d8f9ff1c62f72523f49dab68'

const TWO_PI = 2 * Math.PI

const MANDELBOX = false
const BLOOM = false
const BLOOM_WET = 4.0
const BLOOM_PASSES = 10
const BLOOM_MIN_BRIGHTNESS = 1

// Initialize shell
export default class App {
  constructor (options = {}) {
    let canvas = document.createElement('canvas')
    document.body.appendChild(canvas)
    canvas.style.display = 'none'

    let gl = makeContext(canvas)

    // enable extensions
    var ext = gl.getExtension('OES_standard_derivatives')
    if (!ext) {
      throw new Error('derivatives not supported')
    }

    this.LOOKAT = true

    const preset = {
      offset: {
        x: 1,
        y: 0.78,
        z: 0.669
      },
      d: 5,
      scale: 1.91,
      rot2angle: [0.08, 0, 0],
      cameraAngles: [0, 0, 0]
    }

    this.d = preset.d
    this.cameraRo = vec3.fromValues(0, 0.05, 5)

    // Ray Marching Parameters
    this.epsilon = preset.epsilon || 0.0001

    // Fractal parameters
    this.offset = (preset.offset)
      ? vec3.fromValues(preset.offset.x, preset.offset.y, preset.offset.z)
      : vec3.fromValues(0, 0, 0)
    this.scale = preset.scale
    this.rot2angle = preset.rot2angle || [0, 0, 0]
    this.cameraAngles = preset.cameraAngles || [0, 0, 0]

    // this.setupAnimation(preset)

    this.glInit(gl)

    // Audio
    this.audioFFT = 512
    this.audioTexArray = new Uint8Array(1 * this.audioFFT)
    this.audioNday = ndarray(this.audioTexArray, [this.audioFFT, 1])
    this.audioTex = createTexture(gl, this.audioNday)
    this.pulseGoal = 0

    // Capturing state
    this.capturing = defined(options.capturing, false)

    let tMatCapImg = new Image()
    tMatCapImg.src = './env.jpg'

    let tMatCapImgLoaded = new Promise((resolve, reject) => {
      tMatCapImg.onload = () => {
        this.tMatCap = createTexture(gl, tMatCapImg)
        resolve()
      }
    })

    // this.audioReady = this.setupAudio()
    this.loaded = Promise.all([tMatCapImgLoaded])

    // Scene Rendering
    this.sceneRender = defined(options.sceneRender, this.defaultSceneRender)

    Object.assign(this, {
      canvas,
      gl
    })
  }

  getDimensions () {
    return [dpr * window.innerWidth, dpr * window.innerHeight]
  }

  setupFBOs (gl) {
    let dim = this.getDimensions()
    this.state = [
      createFBO(gl, dim, { depth: false }),
      createFBO(gl, dim, { depth: false }),
      createFBO(gl, dim, { depth: false }) ]

    this.state[0].color.magFilter = gl.LINEAR
    this.state[0].color.minFilter = gl.LINEAR
    this.state[1].color.magFilter = gl.LINEAR
    this.state[1].color.minFilter = gl.LINEAR
    this.state[2].color.magFilter = gl.LINEAR
    this.state[2].color.minFilter = gl.LINEAR
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

    let cameraPosTween = new TWEEN.Tween(ob)
    cameraPosTween
      .to({ y: 0.748, z: 0.934 }, 5 * 1000)
      .onUpdate(updatePos)

    let cameraPosTween2 = new TWEEN.Tween({ x: 0, y: 0.748, z: 0.934 })
    cameraPosTween2
      .to({ y: 0.576, z: 0.752 }, 5 * 1000)
      .easing(TWEEN.Easing.Quadratic.InOut)
      .onUpdate(updatePos)

    let cameraPosTween3 = new TWEEN.Tween({ x: 0.326, y: 0.775, z: 0.752 })
    cameraPosTween3
      .to({ x: 1.016, z: 0.55 }, 5 * 1000)
      .easing(TWEEN.Easing.Quadratic.InOut)
      .onUpdate(updatePos)

    cameraPosTween.chain(cameraPosTween2)
    cameraPosTween2.chain(cameraPosTween3)
    cameraPosTween.start(0)

    // Camera rotation
    function updateRot () {
      self.cameraAngles[0] = this[0]
      self.cameraAngles[1] = this[1]
      self.cameraAngles[2] = this[2]
    }

    let rotObj = [...this.cameraAngles]
    let camRotTween1 = new TWEEN.Tween(rotObj)
    camRotTween1
      .to([0, 0, 0], 0.01 * 1000)
      .onUpdate(updateRot)
      .delay(5 * 1000)
      .easing(TWEEN.Easing.Quadratic.InOut)

    // camRotTween1.start(0)

    // Animation Fractal
    let rotTween1 = new TWEEN.Tween(this.rot2angle)
    rotTween1
      .to([0.777, 0.0, 0], 5 * 1000)
      .easing(TWEEN.Easing.Quadratic.InOut)

    let rotTween2 = new TWEEN.Tween(this.rot2angle)
    rotTween2
      .delay(5 * 1000)
      .to([0.226, 0.0, 5.849 - TWO_PI], 5 * 1000)
      .easing(TWEEN.Easing.Quadratic.InOut)

    rotTween1.chain(rotTween2)
    rotTween1.start(0)

    // Scale Tween
    let scaleTween1 = new TWEEN.Tween(this)
    scaleTween1
      .to({ scale: 5.91 }, 10 * 1000)
      .easing(TWEEN.Easing.Quadratic.InOut)

    // Offset Tween
    let offsetTween1 = new TWEEN.Tween(this.offset)
    offsetTween1
      .delay(10 * 1000)
      .to([
        1,
        -0.069,
        -0.026
      ], 5 * 1000)
      .easing(TWEEN.Easing.Quadratic.InOut)

    offsetTween1.start(0)
  }

  setupAudio () {
    return new Promise((resolve, reject) => {
      SoundCloud({
        client_id: CLIENT_ID,
        song: 'https://soundcloud.com/xlr8r/download-skudge-traveller?in=xlr8r/sets/xlr8rs-top-10-downloads-of-21',
        dark: false,
        getFonts: true
      }, (err, src, data, div) => {
        if (err) {
          reject(err)
          throw err
        }

        // Play the song on
        // a modern browser
        let audio = new Audio()
        audio.crossOrigin = 'Anonymous'
        audio.src = src
        audio.addEventListener('canplay', () => {
          resolve()

          this.audioCtx = new AudioContext()
          let media = this.audioCtx.createMediaElementSource(audio)
          this.analyser = this.audioCtx.createAnalyser()
          // this.analyser.smoothingTimeConstant = 0.5
          this.analyser.fftSize = this.audioFFT
          media.connect(this.analyser)
          this.analyser.connect(this.audioCtx.destination)
          audio.play()
        })

        this.audio = audio

        // Metadata related to the song
        // retrieved by the API.
        console.log(data)
      })
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

  glInit (gl) {
    // Turn off depth test
    gl.disable(gl.DEPTH_TEST)

    // Create fragment shader
    this.shader = createShader(gl, glslify('./vert.glsl'), glslify('./frag.glsl'))
    this.bright = createShader(gl, glslify('./vert.glsl'), glslify('./bright.glsl'))
    this.bloom = createShader(gl, glslify('./vert.glsl'), glslify('./bloom.glsl'))
    this.finalPass = createShader(gl, glslify('./vert.glsl'), glslify('./final-pass.glsl'))

    this.setupFBOs(gl)

    this.shader.attributes.position.location = 0
  }

  kifsM (t = 0) {
    let { offset, scale } = this
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
    fit(canvas, window, dpr)
    let dim = this.getDimensions()

    this.state[0].shape = dim
    this.state[1].shape = dim
    this.state[2].shape = dim
  }

  tick (t) {
    let gl = this.gl

    let dim = this.getDimensions()

    this.shader.bind()

    gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT)
    gl.clearColor(0, 0, 0, 1)
    gl.viewport(0, 0, dim[0], dim[1])

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
    TWEEN.update(t)

    if (this.tMatCap) {
      this.shader.uniforms.tMatCap = this.tMatCap.bind(0)
    }

    // Update audio
    if (this.analyser) {
      this.analyser.getByteTimeDomainData(this.audioTexArray)
      this.audioTex.setPixels(this.audioNday)
      this.shader.uniforms.audioTexture = this.audioTex.bind(1)
    }

    this.shader.uniforms.epsilon = this.epsilon

    let updates = this.getCamera(t)
    this.shader.uniforms.cameraRo = updates[0]
    this.shader.uniforms.cameraMatrix = (updates[1])

    this.shader.uniforms.kifsM = this.kifsM(t)
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
      drawTriangle(gl)

      // Vertical Blur
      let prev = this.state[2].color[0]
      this.state[1].bind()

      this.bloom.uniforms.buffer = prev.bind(2)
      this.bloom.uniforms.resolution = dim
      this.bloom.uniforms.direction = [0, 1]
      drawTriangle(gl)
    }

    // Additive blending
    gl.bindFramebuffer(gl.FRAMEBUFFER, null)
    this.finalPass.bind()
    this.finalPass.uniforms.base = this.state[0].color[0].bind(3)
    this.finalPass.uniforms.buffer = this.state[1].color[0].bind(4)
    this.finalPass.uniforms.resolution = dim
    this.finalPass.uniforms.time = this.getTime(t)
    this.finalPass.uniforms.wet = BLOOM_WET
    drawTriangle(gl)
  }

  getTime (t) {
    return window.time || t / 1000
  }

  defaultSceneRender (_, t) {
    drawTriangle(this.gl)
  }

  render (t) {
    let { shader, gl } = this

    if (BLOOM) {
      this.state[0].bind()
    }

    shader.uniforms.time = this.getTime(t)
    shader.uniforms.BLOOM = BLOOM
    this.sceneRender(shader, t)

    if (BLOOM) {
      this.bloomBlur(gl, t)
    }
  }

  run () {
    this.canvas.style.display = null
    this.resize()

    this.enableEvents()
  }

  stop () {
    this.canvas.style.display = 'none'
    this.disposeEvents()
  }
}
