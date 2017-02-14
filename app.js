const glslify = require('glslify')

import createShader from 'gl-shader'
import createTexture from 'gl-texture2d'

import ShaderVREffect from 'shader-vr-effect'
import ShaderVROrbitControls from 'shader-vr-orbit-controls'
import WebVRManager from 'shader-webvr-manager'

import fit from 'canvas-fit'
import TWEEN from 'tween.js'
import makeContext from 'gl-context'
import { isAndroid, rot4 } from './utils'
import CCapture from 'ccapture.js'

import assign from 'object-assign'
import defined from 'defined'
import { vec3, mat4 } from 'gl-matrix'
import presets from './presets.json'

const dpr = Math.min(2, defined(window.devicePixelRatio, 1))

const fr = 60
let captureTime = 0
const secondsLong = 38

const capturing = false
const LOOKAT = false

let capturer = {}
if (capturing) {
  capturer = new CCapture({
    format: 'jpg',
    framerate: fr,
    name: 'kifs-icosa-spelunking',
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
        x: 1,
        y: 0,
        z: PHI
      },
      d: 3,
      scale: 2,
      rot2angle: 0
    }

    // Ray Marching Parameters
    this.epsilon = preset.epsilon || 0.005

    // Fractal parameters
    this.offset = vec3.fromValues(preset.offset.x, preset.offset.y, preset.offset.z)
    this.d = preset.d
    this.scale = preset.scale
    this.rot2angle = preset.rot2angle || [0, 0, 0]
    this.cameraAngles = preset.cameraAngles || [0, 0, 0]

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

    eps1.start(0)

    this.cameraRo = vec3.fromValues(0, 0, this.d) // vec3.fromValues(.64, .32, 1.45) // 

    // Camera location animation
    let tween1 = new TWEEN.Tween(this.cameraRo)
    tween1
      .to([.45, .25, 1.6], 7000)
      .easing(TWEEN.Easing.Quadratic.InOut)
    let tween2 = new TWEEN.Tween(this.cameraRo)
    tween2
      .delay(1000)
      .to([.62, .31, 1.47], 5500)
      .easing(TWEEN.Easing.Quadratic.InOut)
    let tween3 = new TWEEN.Tween(this.cameraRo)
    tween3
      .to([.641, .323, 1.446], 5 * 1000)
      .easing(TWEEN.Easing.Quadratic.InOut)
    let tween4 = new TWEEN.Tween(this.cameraRo)
    tween4
      .delay(1000)
      .to([0, 0, this.d], 11 * 1000)
      .easing(TWEEN.Easing.Quadratic.InOut)

    tween1.chain(tween2)
    tween2.chain(tween3)
    tween3.chain(tween4)

    tween1.start(0)

    // Camera rotation animation
    let camRotTween1 = new TWEEN.Tween(this.cameraAngles)
    camRotTween1
      .delay(7000)
      .to([Math.PI / 7, -Math.PI / 8, 0], 2000)
      .easing(TWEEN.Easing.Quadratic.InOut)
    let camRotTween2 = new TWEEN.Tween(this.cameraAngles)
    camRotTween2
      .delay(2000)
      .to([-2 * Math.PI + 5.535, - 2 * Math.PI + 5.806, -2 * Math.PI + 5.041], 5500)
      .easing(TWEEN.Easing.Quadratic.InOut)
    let camRotTween3 = new TWEEN.Tween(this.cameraAngles)
    camRotTween3
      .delay(1000)
      .to([0, 0, 0], 9 * 1000)
      .easing(TWEEN.Easing.Quadratic.InOut)

    camRotTween1.chain(camRotTween2)
    camRotTween2.chain(camRotTween3)

    camRotTween1.start(0)

    // Animation Fractal
    let rotTween1 = new TWEEN.Tween(this.rot2angle)
    rotTween1
      .to([-1 / 4 * Math.PI, 0, -1 / 4 * Math.PI], 10000)
      .easing(TWEEN.Easing.Quadratic.InOut)
      .delay(5000)
    let rotTween2 = new TWEEN.Tween(this.rot2angle)
    rotTween2
      .to([1 / 3 * Math.PI, 0, 1 / 3 * Math.PI], 10000)
      .easing(TWEEN.Easing.Cubic.InOut)
    let rotTween3 = new TWEEN.Tween(this.rot2angle)
    rotTween3
      .to([0, 0, 0], 5000)
      .easing(TWEEN.Easing.Cubic.InOut)

    rotTween1.chain(rotTween2)
    rotTween2.chain(rotTween3)

    // rotTween1.start(0)

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
      running: false,
    })

    this.setupStage()
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

    let img = document.createElement('img')
    img.src = 'env.jpg'

    // Create fragment shader
    this.shader = createShader(gl, glslify('./vert.glsl'), glslify('./frag.glsl'))

    this.shader.bind()

    img.onload = () => {
      this.texture = createTexture(gl, img)
    }

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
    navigator.getVRDisplays().then((displays) => {
      if (displays.length > 0) {
        this.vrDisplay = displays[0]
        this.effect.setVRDisplay(this.vrDisplay)
        this.controls.setVRDisplay(this.vrDisplay)

        if (this.running && !this.currentRAF) {
          this.currentRAF = this.vrDisplay.requestAnimationFrame(this.tick.bind(this))
        }
      }
    })
  }

  resize (e) {
    let { effect, canvas } = this
    let scale = 1
    fit(canvas, window, dpr * scale)

    effect.setSize(dpr * scale * window.innerWidth, dpr * scale * window.innerHeight)
  }

  tick (t) {
    t = capturing ? currentTime + 1000 / fr : t
    currentTime = t

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

  render (t) {
    let { shader, manager, controls } = this

    shader.uniforms.time = t / 1000
    shader.uniforms.texture = this.texture

    manager.render(shader, t)
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

    if (this.vrDisplay && !this.currentRAF) {
      this.currentRAF = this.vrDisplay.requestAnimationFrame(this.tick.bind(this))
    }
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
