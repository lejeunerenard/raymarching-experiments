const glslify = require('glslify')

import createShader from 'gl-shader'
import createTexture from 'gl-texture2d'

import ShaderVREffect from 'shader-vr-effect'
import ShaderVROrbitControls from 'shader-vr-orbit-controls'
import WebVRManager from 'shader-webvr-manager'

import fit from 'canvas-fit'
import makeContext from 'gl-context'
import { isAndroid } from './utils'
import CCapture from 'ccapture.js'

import assign from 'object-assign'
import defined from 'defined'

const dpr = Math.min(2, defined(window.devicePixelRatio, 1))

const fr = 120
let captureTime = 0
const secondsLong = 20

const capturing = false

let capturer = {}
if (capturing) {
  capturer = new CCapture({
    format: 'jpg',
    framerate: fr,
    name: 'kelido-crystal',
    autoSaveTime: 10,
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

    this.setupStage()
  }

  enableEvents () {
    return
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
    img.onload = () => {
      this.texture = createTexture(gl, img)
    }

    this.shader.attributes.position.location = 0
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
    t = currentTime + 1000 / fr
    currentTime = t

    this.update(t)
    this.render(t)

    this.currentRAF = this.vrDisplay.requestAnimationFrame(this.tick.bind(this))
  }

  update () {
  }

  render (t) {
    let { shader, manager, controls } = this

    shader.bind()

    shader.uniforms.time = t / 1000 / 4
    shader.uniforms.texture = this.texture
    controls.update(shader)

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
