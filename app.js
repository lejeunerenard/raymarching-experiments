import createShader from 'gl-shader'
import createTexture from 'gl-texture2d'
import createFBO from 'gl-fbo'

import ndarray from 'ndarray'
import TWEEN from 'tween.js'
import fill from 'ndarray-fill'
import makeContext from 'gl-context'
import { rot4 } from './utils'
import drawTriangle from 'a-big-triangle'

import defined from 'defined'
import assert from 'assert'
import { vec3, mat4 } from 'gl-matrix'
import { getLuminance, getColorWFixedLuminance } from './luminance'

// import convert from './svg-to-glsl.js'

const dpr = 1.0 / window.devicePixelRatio

const TWO_PI = 2 * Math.PI
// const PHI = (1 + Math.sqrt(5)) / 2

const MANDELBOX = false
const BLOOM = false
const BLOOM_PASSES = 4
const BLOOM_RADIUS = BLOOM_PASSES - 1
const BLOOM_WET = 2.0
const BLOOM_MIN_BRIGHTNESS = 0.85

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
    this.SHOW_SVG_SDF = false

    this.presets = {}
    const preset = {
      offset: {
        x: 0.492,
        y: 0.218,
        z: 0.197
      },
      d: 0.52,
      scale: 1.3325,
      rot2angle: [0.076, 0.716, 0.451],
      cameraAngles: [0.05, 0.085, 0.289]
    }

    this.d = preset.d
    this.cameraRo = vec3.fromValues(-1.65, 1.65, 1.65)
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

    this.angle1C = -0.2654
    this.angle2C = 0.449
    this.angle3C = 0.89

    // this.setupAnimation(preset)

    this.glInit(gl)

    // Audio
    this.audioFFT = 32
    this.audioTexArray = new Uint8Array(1 * this.audioFFT)
    this.audioNday = ndarray(this.audioTexArray, [this.audioFFT, 1])
    this.audioTex = createTexture(gl, this.audioNday)
    this.pulseGoal = 0

    // -- Shapes --
    this.shapeOptions = {
      cube: 0,
      sphere: 1,
      cigar: 2,
      pyramid: 3,
      torus: 4,
      octahedron: 5,
      none: -1
    }
    this.shapeMode = 'cube'
    this.shapeScale = [1, 1, 1]
    this._shape2DSDFTextures = {}

    // Capturing state
    this.capturing = defined(options.capturing, false)

    this.loaded = Promise.resolve()
      // .then(() => {
      //   const { glsl, metadata } = convert(require('./7.svg.js'))
      //   const fbo = this.generateSVGTexture(glsl, metadata.viewBox)
      //   this.add2DSDFTexture('year-7', fbo.color[0])
      // })
      .then(() => {
        this.setupAudio()
      })

    // Scene Rendering
    this.sceneRenderer = options.sceneRenderer

    this.totalTime = 0

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

  renderSVGSDF (svgToSDF, viewBox = { x1: 0, y1: 0, x2: 100, y2: 100 }, fbo, gl = this.gl) {
    const t = 0
    const MAX_TEXTURE_SIZE = gl.getParameter(gl.MAX_TEXTURE_SIZE)
    let dim = [MAX_TEXTURE_SIZE / 2, MAX_TEXTURE_SIZE / 2]

    let svgSDF
    // By default create a FBO to return
    if (!fbo) {
      try {
        svgSDF = createFBO(gl, dim, { preferFloat: true, depth: false })
      } catch (e) {
        svgSDF = createFBO(gl, dim, { depth: false })
      }
    } else {
      svgSDF = fbo
    }
    svgSDF.color[0].wrap = [gl.CLAMP_TO_EDGE, gl.CLAMP_TO_EDGE]
    svgSDF.color[0].magFilter = gl.LINEAR
    svgSDF.color[0].minFilter = gl.LINEAR
    svgSDF.color[0].mipSamples = 12

    svgToSDF.bind()
    svgToSDF.uniforms.resolution = dim
    svgToSDF.uniforms.scale = Math.min(dim[0], dim[1]) / Math.max(viewBox.x2 - viewBox.x1, viewBox.y2 - viewBox.y1)
    svgToSDF.uniforms.epsilon = this.epsilon

    let updates = this.getCamera(t)
    svgToSDF.uniforms.cameraRo = this.SHOW_SVG_SDF ? updates[0] : [0, 0, 1.25]
    svgToSDF.uniforms.cameraMatrix = (updates[1])

    if (this.SHOW_SVG_SDF) {
      gl.bindFramebuffer(gl.FRAMEBUFFER, null)
    } else {
      svgSDF.bind()
    }

    drawTriangle(gl)

    return svgSDF
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
      .to({ x: self.cameraRo[0], y: -0.1, z: self.cameraRo[2] }, 5 * 1000)
      .easing(TWEEN.Easing.Quadratic.Out)
      .onUpdate(updatePos)

    let cameraPosTween2 = new TWEEN.Tween(ob)
    cameraPosTween2
      .to({ x: self.cameraRo[0], y: self.cameraRo[1], z: self.cameraRo[2] }, 5 * 1000)
      .easing(TWEEN.Easing.Quadratic.Out)
      .onUpdate(updatePos)

    cameraPosTween.chain(cameraPosTween2)
    cameraPosTween2.chain(cameraPosTween)
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
      .to([2.144, this.rot2angle[1], this.rot2angle[2]], 6 * 1000)
      .easing(TWEEN.Easing.Quadratic.InOut)
    // let rotTween2 = new TWEEN.Tween(this.rot2angle)
    // rotTween2
    //   .delay(1 * 1000)
    //   .to([0.125, this.rot2angle[1], 0.985], 14 * 1000)
    //   .easing(TWEEN.Easing.Quadratic.InOut)
    let rotTween3 = new TWEEN.Tween(this.rot2angle)
    rotTween3
      .delay(2 * 1000)
      .to([...this.rot2angle], 6 * 1000)
      .easing(TWEEN.Easing.Quadratic.InOut)

    rotTween1.chain(rotTween3)
    // rotTween2.chain(rotTween3)
    rotTween3.chain(rotTween1)

    // rotTween1.start(0)

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
      .to([ this.offset[0], 0.539, this.offset[2] ], 10 * 1000)
      .easing(TWEEN.Easing.Quadratic.InOut)
    let offsetTween2 = new TWEEN.Tween(this.offset)
    offsetTween2
      .to([ this.offset[0], this.offset[1], this.offset[2] ], 10 * 1000)
      .easing(TWEEN.Easing.Quadratic.InOut)

    offsetTween1.chain(offsetTween2)
    offsetTween2.chain(offsetTween1)

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

    this.generateSVGTexture('', undefined, gl)
    gl.bindFramebuffer(gl.FRAMEBUFFER, null) // Undo whatever generateSVGTexture bound

    let blankSDF = ndarray([], [2, 2])
    // Initialize at max distance
    fill(blankSDF, (i) => 1.0)
    this.defaultSDFTexture = createTexture(gl, blankSDF)

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

  generateSVGTexture (svgPaths, viewBox = { x1: 0, y1: 0, x2: 100, y2: 100 }, gl = this.gl) {
    const svgShader = require('./shaders/svgSDF')(svgPaths)
    this.svgToSDF = createShader(gl, svgShader.vertex, svgShader.fragment)
    this.svgViewBox = viewBox
    return this.renderSVGSDF(this.svgToSDF, this.svgViewBox, null, gl)
  }

  add2DSDFTexture (name, texture) {
    this.shapeOptions[name] = name
    this._shape2DSDFTextures[name] = {
      tex: texture,
      filename: name,
      _lastFilename: name,
      isVideo: false,
      asset: null
    }
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
    if (this.SHOW_SVG_SDF) {
      if (this.debugSVGFBO) {
        this.renderSVGSDF(this.svgToSDF, this.svgViewBox, this.debugSVGFBO)
      } else {
        this.debugSVGFBO = this.renderSVGSDF(this.svgToSDF, this.svgViewBox)
      }
    }
  }

  getCamera (t) {
    t = this.getTime(t)
    let cameraMatrix = mat4.create()

    const origin = vec3.fromValues(0, 0, 0)
    const unitX = vec3.fromValues(1, 0, 0)
    const unitY = vec3.fromValues(0, 1, 0)
    const unitZ = vec3.fromValues(0, 0, 1)

    // Camera Rotation
    const cameraRo = vec3.clone(this.cameraRo)
    // vec3.rotateY(cameraRo, cameraRo, origin, TWO_PI * t / this.totalTime)

    // LookAt
    if (this.LOOKAT) {
      mat4.lookAt(cameraMatrix, cameraRo, origin, unitY)
    } else {
      const angleX = this.cameraAngles[0]
      const axisX = unitX
      mat4.multiply(cameraMatrix, rot4(axisX, angleX), cameraMatrix)

      // Y-centric
      const angleY = this.cameraAngles[1]
      const axisY = unitY
      mat4.multiply(cameraMatrix, rot4(axisY, angleY), cameraMatrix)

      // Z-centric
      const angleZ = this.cameraAngles[2]
      const axisZ = unitZ
      mat4.multiply(cameraMatrix, rot4(axisZ, angleZ), cameraMatrix)
    }

    this.cameraMatrix = cameraMatrix
    return [cameraRo, cameraMatrix]
  }

  update (t) {
    t = (window.time !== undefined) ? window.time : t

    TWEEN.update(t)

    this.shader.uniforms.epsilon = this.epsilon

    // Shape
    const SDF_TEX_LOC = 2
    const isSVG = true
    const svgName = 'year-7'
    let sdfTexture = this._shape2DSDFTextures[svgName]
    if (isSVG && sdfTexture) {
      this.shader.uniforms.sdf2DTexture = sdfTexture.tex.bind(SDF_TEX_LOC)
    } else {
      this.shader.uniforms.sdf2DTexture = this.defaultSDFTexture.bind(SDF_TEX_LOC)
    }

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
      this.bloom.uniforms.direction = [BLOOM_RADIUS - i, 0]
      this.bloom.uniforms.time = this.getTime(t)
      drawTriangle(gl)

      // Vertical Blur
      let prev = this.state[2].color[0]
      this.state[1].bind()

      this.bloom.uniforms.buffer = prev.bind(2)
      this.bloom.uniforms.resolution = dim
      this.bloom.uniforms.direction = [0, BLOOM_RADIUS - i]
      this.bloom.uniforms.time = this.getTime(t)
      drawTriangle(gl)
    }
  }

  finalPassRender (gl, t) {
    let dim = this.getDimensions()

    let base = this.state[0].color[0]

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

    // Always render to FBO as finalPassRender will render to final frame buffer
    this.state[0].bind()

    shader.uniforms.time = this.getTime(t)
    shader.uniforms.BLOOM = BLOOM
    this.sceneRenderer(shader, t)

    if (BLOOM) {
      this.bloomBlur(gl, t)
    }

    this.finalPassRender(gl, t)
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
