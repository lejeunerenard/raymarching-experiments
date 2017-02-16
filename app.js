const glslify = require('glslify')

import createShader from 'gl-shader'
import createTexture from 'gl-texture2d'
import createFBO from 'gl-fbo'

import ShaderVREffect from 'shader-vr-effect'
import ShaderVROrbitControls from 'shader-vr-orbit-controls'
import WebVRManager from 'shader-webvr-manager'

import fit from 'canvas-fit'
import TWEEN from 'tween.js'
import makeContext from 'gl-context'
import { isAndroid, rot4 } from './utils'
import { cameraOrbit } from './camera-tweens'
import CCapture from 'ccapture.js'
import SoundCloud from 'soundcloud-badge'
import Analyser from 'gl-audio-analyser'
import drawTriangle from 'a-big-triangle'

import assign from 'object-assign'
import defined from 'defined'
import { vec3, mat4 } from 'gl-matrix'
import presets from './presets.json'

const dpr = Math.min(2, defined(window.devicePixelRatio, 1))
const CLIENT_ID = 'ded451c6d8f9ff1c62f72523f49dab68'

const fr = 25
let captureTime = 0
const secondsLong = 10

const capturing = false
const LOOKAT = true

let capturer = {}
if (capturing) {
  capturer = new CCapture({
    format: 'jpg',
    framerate: fr,
    name: 'kifs-swanky-heart',
    autoSaveTime: 5,
    startTime: captureTime,
    timeLimit: secondsLong,
    verbose: true
  })
}

let currentTime = captureTime * 1000
window.capturer = capturer
let winSetTimeout = window.setTimeout
let winClearTimeout = window.clearTimeout
let winSetInterval = window.setInterval
let winclearInterval = window.clearInterval
let winRequestAnimationFrame = window.requestAnimationFrame
let winProfNow = window.performance.now

const PHI = (1+Math.sqrt(5))/2;

// Initialize shell
export default class App {
  constructor (options = {}) {
    let canvas = document.createElement('canvas')
    document.body.appendChild(canvas)
    canvas.style.display = 'none'
    if (!isAndroid()) {
      canvas.addEventListener('touchstart', function(e) { e.preventDefault() })
    }

    let gl = makeContext(canvas)

    const preset = {
      offset: {
        x: 1.883,
        y: .669,
        z: 1.763
      },
      d: 4.5,
      scale: 2.56,
      rot2angle: 0
    }

    this.d = preset.d
    this.cameraRo = vec3.fromValues(0, 0, this.d)

    // Ray Marching Parameters
    this.epsilon = preset.epsilon || 0.003

    // Fractal parameters
    this.offset = (preset.offset)
      ? vec3.fromValues(preset.offset.x, preset.offset.y, preset.offset.z)
      : vec3.fromValues(0, 0, 0)
    this.scale = preset.scale
    this.rot2angle = preset.rot2angle || [0, 0, 0]
    this.cameraAngles = preset.cameraAngles || [0, 0, 0]

    this.setupAnimation(preset)

    this.glInit(gl)

    let effect = new ShaderVREffect(gl)
    let controls = new ShaderVROrbitControls(gl)

    let params = {
      hideButton: true,
      isUndistorted: false
    }
    let manager = new WebVRManager({ domElement: canvas }, effect, params)
    let vrDisplay = undefined

    assign(this, {
      canvas,
      gl,
      effect,
      controls,
      manager,
      vrDisplay,
      currentRAF: null,
      running: false
    })

    this.stageReady = this.setupStage()
    this.loaded = Promise.all([this.stageReady])
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
    // Epsilon Animation
    let eps1 = new TWEEN.Tween(this)
    eps1
      .to({ epsilon: 0.0001 }, 7000)
      .easing(TWEEN.Easing.Quadratic.InOut)
    let eps2 = new TWEEN.Tween(this)
    eps2
      .to({ epsilon: 0.00006 }, 6000)
      .easing(TWEEN.Easing.Quadratic.InOut)
    let eps3 = new TWEEN.Tween(this)
    eps3
      .to({ epsilon: 0.000009 }, 5000)
      .easing(TWEEN.Easing.Quadratic.InOut)
    let eps4 = new TWEEN.Tween(this)
    eps4
      .delay(2000)
      .to({ epsilon: 0.005 }, 10000)
      .easing(TWEEN.Easing.Quadratic.InOut)

    eps1.chain(eps2)
    eps2.chain(eps3)
    eps3.chain(eps4)

    // eps1.start(0)

    this.cameraRo = vec3.fromValues(0, 0, this.d)

    // Camera location animation
    let posRot = [0, 0]

    let cameraPosTween =
      cameraOrbit(this.cameraRo, this.d, posRot, [Math.PI / 8, Math.PI / 4], 7 * 1000)
    let cameraPosTween2 =
      cameraOrbit(this.cameraRo, this.d, [-Math.PI / 3, Math.PI / 2], [-Math.PI / 3, 0], 7 * 1000)
    let cameraPosTween3 =
      cameraOrbit(this.cameraRo, this.d, [-Math.PI / 4, -Math.PI * 3 / 4], [Math.PI / 4, Math.PI / 2], 7 * 1000)

    cameraPosTween.chain(cameraPosTween2)
    cameraPosTween2.chain(cameraPosTween3)
    cameraPosTween3.chain(cameraPosTween)
    cameraPosTween.start(0)

    // Animation Fractal
    let rotTween1 = new TWEEN.Tween(this.rot2angle)
    rotTween1
      .to([-1 / 3 * Math.PI, .5, -1 / 4 * Math.PI], 10000)
      .easing(TWEEN.Easing.Quadratic.InOut)
      .delay(1000)
    let rotTween2 = new TWEEN.Tween(this.rot2angle)
    rotTween2
      .to([0, 0, 0], 10000)
      .easing(TWEEN.Easing.Cubic.InOut)

    rotTween1.chain(rotTween2)
    rotTween2.chain(rotTween1)

    rotTween1.start(0)
  }

  setupAudio () {
    return new Promise((resolve, reject) => {
      SoundCloud({
        client_id: CLIENT_ID,
        song: 'https://soundcloud.com/shang-lin/swanky',
        dark: true,
        getFonts: true
      }, (err, src, data, div)  => {
        if (err) {
          reject(err)
          throw err
        }

        // Play the song on
        // a modern browser
        let audio = new Audio
        audio.crossOrigin = 'Anonymous'
        audio.src = src
        audio.addEventListener('canplay', () => {
          console.log('playing!')
          this.analyser = Analyser(this.gl, audio)
          audio.play()
        })

        this.audio = audio

        resolve()

        // Metadata related to the song
        // retrieved by the API.
        console.log(data)
      })
    })
  }

  enableEvents () {
    if (capturing) return
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
    let _kifsM = mat4.fromValues(
      scale, 0,     0,     -offset[0] * (scale - 1),
      0,     scale, 0,     -offset[1] * (scale - 1),
      0,     0,     scale, -offset[2] * (scale - 1),
      0,     0,     0,     1)

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

  // Get the HMD, and if we're dealing with something that specifies
  // stageParameters, rearrange the scene.
  setupStage () {
    return navigator.getVRDisplays().then((displays) => {
      if (displays.length > 0) {
        this.vrDisplay = displays[0]
        this.effect.setVRDisplay(this.vrDisplay)
        this.controls.setVRDisplay(this.vrDisplay)
      }
    })
  }

  resize (e) {
    let { effect, canvas } = this
    let scale = 1
    fit(canvas, window, dpr * scale)
    let dim = this.getDimensions()

    effect.setSize(scale * dim[0], scale * dim[1])
    this.state[0].shape = dim
    this.state[1].shape = dim
    this.state[2].shape = dim
  }

  tick (t) {
    let gl = this.gl

    let dim = this.getDimensions()

    t = capturing ? currentTime + 1000 / fr : t
    currentTime = t

    this.shader.bind()

    gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT)
    gl.clearColor(0, 0, 0, 1)
    gl.viewport(0, 0, dim[0], dim[1])

    this.update(t)
    this.render(t)

    this.currentRAF = this.vrDisplay.requestAnimationFrame(this.tick.bind(this))
  }

  getCamera (t) {
    t /= 1000
    let cameraMatrix = mat4.create() 

    // LookAt
    if (LOOKAT) {
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

    this.shader.uniforms.epsilon = this.epsilon

    this.controls.update(this.shader)

    let updates = this.getCamera(t)
    this.shader.uniforms.cameraRo = updates[0]
    this.shader.uniforms.cameraMatrix = (updates[1])

    this.shader.uniforms.kifsM = this.kifsM(t)
  }

  blur (gl) {
    let dim = this.getDimensions()

    // Brightness pass
    let base = this.state[0].color[0]
    this.state[1].bind()
    this.bright.bind()
    this.bright.uniforms.minBright = .2
    this.bright.uniforms.buffer = base.bind(0)
    this.bright.uniforms.resolution = dim
    drawTriangle(gl)

    for (let i = 0; i < 5; i++) {
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

      //this.bloom.bind()
      this.bloom.uniforms.buffer = prev.bind(2)
      this.bloom.uniforms.resolution = dim
      this.bloom.uniforms.direction = [0, 1]
      drawTriangle(gl)
    }

    // Additive blending
    gl.bindFramebuffer(gl.FRAMEBUFFER, null)
    this.finalPass.bind()
    this.finalPass.uniforms.base = base.bind(0)
    this.finalPass.uniforms.buffer = this.state[1].color[0].bind(1)
    this.finalPass.uniforms.resolution = dim
    drawTriangle(gl)
  }

  render (t) {
    let { shader, manager, controls, gl } = this

    this.state[0].bind()
    shader.uniforms.time = t / 1000
    manager.render(shader, t)

    this.blur(gl)

    capturing && capturer.capture(this.canvas)
  }

  run () {
    this.canvas.style.display = null
    this.resize()

    if (this.manager.mode !== WebVRManager.Modes.VR) {
      this.manager.button.setVisibility(true)
    }

    this.enableEvents()
    this.manager.enableEvents()

    this.running = true
    capturing && capturer.start()
    window.setTimeout = winSetTimeout
    window.clearTimeout = winClearTimeout
    window.setInterval = winSetInterval
    window.clearInterval = winclearInterval
    window.requestAnimationFrame = winRequestAnimationFrame
    window.performance.now = winProfNow

    this.loaded.then(() => {
      if (this.vrDisplay && this.running && !this.currentRAF) {
        this.currentRAF = this.vrDisplay.requestAnimationFrame(this.tick.bind(this))
      }
    })
  }

  stop () {
    this.canvas.style.display = 'none'
    this.manager.button.setVisibility(false)
    if (this.currentRAF) {
      this.vrDisplay.cancelAnimationFrame(this.currentRAF)
      this.currentRAF = null
    }
    this.running = false
    this.disposeEvents()
    this.manager.disposeEvents()
  }
}
