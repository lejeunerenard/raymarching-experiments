import WebVRManager from 'shader-webvr-manager'
import ShaderVREffect from 'shader-vr-effect'
import ShaderVROrbitControls from 'shader-vr-orbit-controls'

import App from './app'
import { name } from './info.json'

const fr = 60
const captureTime = 0 * 5
const secondsLong = 6
const capturing = false

const FOV = 70

let app = new App()
window.app = app

app.width = 1080
app.height = 1920

// window.time = 0.3
const still = false

let capturer = {}
if (capturing) {
  let massagedName = name.replace(/ /g, '-')
  massagedName = massagedName.replace(/'/g, '')
  massagedName = massagedName.toLowerCase()

  let filename = massagedName + '-render2'
  console.log('filename', filename)
  capturer = new CCapture({
    format: 'jpg',
    framerate: fr,
    name: filename,
    autoSaveTime: 5,
    quality: 98,
    startTime: captureTime,
    timeLimit: secondsLong,
    verbose: true
  })

  let currentTime = captureTime * 1000
  window.capturer = capturer
  let winSetTimeout = window.setTimeout
  let winClearTimeout = window.clearTimeout
  let winSetInterval = window.setInterval
  let winclearInterval = window.clearInterval
  let winRequestAnimationFrame = window.requestAnimationFrame
  let winProfNow = window.performance.now

  let run = () => {
    capturer.start()
    window.setTimeout = winSetTimeout
    window.clearTimeout = winClearTimeout
    window.setInterval = winSetInterval
    window.clearInterval = winclearInterval
    window.requestAnimationFrame = winRequestAnimationFrame
    window.performance.now = winProfNow
    app.run()
  }

  let vrDisplay, effect, controls, currentRAF, manager

  effect = new ShaderVREffect(app.gl)
  effect.fov = FOV
  controls = new ShaderVROrbitControls(app.gl)

  let params = {
    hideButton: true,
    isUndistorted: false
  }
  manager = new WebVRManager({ domElement: app.canvas }, effect, params)
  app.sceneRender = manager.render.bind(manager)

  let resize = () => {
    let dim = app.getDimensions()
    effect.setSize(dim[0], dim[1])
  }
  // window.addEventListener('resize', resize)
  resize()

  navigator.getVRDisplays().then((displays) => {
    if (displays.length > 0) {
      vrDisplay = displays[0]
      effect.setVRDisplay(vrDisplay)
      controls.setVRDisplay(vrDisplay)
    }

    if (manager.mode !== WebVRManager.Modes.VR) {
      manager.button.setVisibility(true)
    }

    manager.enableEvents()
    run()
  })

  let tick = () => {
    let t = currentTime + 1000 / fr
    currentTime = t

    app.tick(t)
    capturer.capture(app.canvas)

    if (currentTime <= 1000 * (secondsLong + captureTime) + 1000 / fr) {
      window.setTimeout(() => {
        if (!still) {
          currentRAF = vrDisplay.requestAnimationFrame(tick)
        }
      }, 500)
    } else {
      window.setTimeout(() => {
        console.log('done sending message')
        var xmlHTTP = new XMLHttpRequest()
        xmlHTTP.open('GET', 'http://localhost:7321/', false)
        xmlHTTP.send(null)
        console.log('response text', xmlHTTP.responseText)
      }, 250)
    }
  }

  app.loaded.then(() => {
    if (vrDisplay && !currentRAF) {
      currentRAF = vrDisplay.requestAnimationFrame(tick)
    }
  })
} else {
  let gui = new dat.GUI()
  gui.remember(app)
  gui.closed = true

  let offsetF = gui.addFolder('Offset')
  offsetF.add(app.offset, '0', -5, 5).step(0.001).listen()
  offsetF.add(app.offset, '1', -5, 5).step(0.001).listen()
  offsetF.add(app.offset, '2', -5, 5).step(0.001).listen()

  gui.add(app, 'd', 0, 20).step(0.01).listen()
  gui.add(app, 'scale', -6, 6).step(0.01).listen()
  gui.add(app, 'epsilon', 0.0000001, 0.05).step(0.0000001).listen()

  let angleCF = gui.addFolder('Angle Coefficients')
  angleCF.add(app, 'angle1C', -10, 10).step(0.001).listen()
  angleCF.add(app, 'angle2C', 0, 15).step(0.001).listen()
  angleCF.add(app, 'angle3C', -5, 100).step(0.001).listen()

  let rotationF = gui.addFolder('Rotation')
  rotationF.add(app.rot2angle, '0', 0, 2 * Math.PI).step(0.001).listen()
  rotationF.add(app.rot2angle, '1', 0, 2 * Math.PI).step(0.001).listen()
  rotationF.add(app.rot2angle, '2', 0, 2 * Math.PI).step(0.001).listen()

  let cameraF = gui.addFolder('Camera')

  let cameraPosF = cameraF.addFolder('Position')
  cameraPosF.add(app.cameraRo, '0', -20, 20).step(0.01).listen()
  cameraPosF.add(app.cameraRo, '1', -20, 20).step(0.01).listen()
  cameraPosF.add(app.cameraRo, '2', -20, 20).step(0.01).listen()

  let cameraRotF = cameraF.addFolder('Rotation')
  cameraRotF.add(app.cameraAngles, '0', -Math.PI, Math.PI).step(0.001).listen()
  cameraRotF.add(app.cameraAngles, '1', -Math.PI, Math.PI).step(0.001).listen()
  cameraRotF.add(app.cameraAngles, '2', -Math.PI, Math.PI).step(0.001).listen()

  let colorsF = gui.addFolder('Colors')
  colorsF.addColor(app, 'colors1').listen()
  colorsF.addColor(app, 'colors2').listen()

  // ----- Setup -----

  // VR Setup
  let vrDisplay, effect, controls, currentRAF, manager

  effect = new ShaderVREffect(app.gl)
  effect.fov = FOV
  controls = new ShaderVROrbitControls(app.gl)

  let params = {
    hideButton: true,
    isUndistorted: false
  }
  manager = new WebVRManager({ domElement: app.canvas }, effect, params)

  let resize = () => {
    let dim = app.getDimensions()
    effect.setSize(dim[0], dim[1])
  }
  window.addEventListener('resize', resize)
  resize()

  navigator.getVRDisplays().then((displays) => {
    if (displays.length > 0) {
      vrDisplay = displays[0]
      effect.setVRDisplay(vrDisplay)
      controls.setVRDisplay(vrDisplay)
    }

    if (manager.mode !== WebVRManager.Modes.VR) {
      manager.button.setVisibility(true)
    }

    manager.enableEvents()
  })

  // Run
  app.sceneRender = manager.render.bind(manager)
  app.run()

  let tick = (t) => {
    if (!still) {
      currentRAF = vrDisplay.requestAnimationFrame(tick)
    }

    controls.update(app.shader)
    app.tick(t)
  }

  app.loaded.then(() => {
    if (vrDisplay && !currentRAF) {
      currentRAF = vrDisplay.requestAnimationFrame(tick)
    }
  })
}
