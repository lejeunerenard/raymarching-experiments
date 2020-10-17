import DefaultSceneRenderer from './default-scene-renderer'
import App from './app'

import { name } from './info.json'

const fr = 60
const captureTime = 0 * 5
const secondsLong = 30
const capturing = false

let app = new App()
window.app = app

app.width = 1080
app.height = 1920

let sceneRenderer = new DefaultSceneRenderer(app.gl)
app.sceneRenderer = sceneRenderer.render.bind(sceneRenderer)

// Setup resizing
let resize = () => {
  let dim = app.getDimensions()
  sceneRenderer.resize(dim[0], dim[1])
}
resize()

// window.time = 0.3
const still = false

if (capturing) {
  let capturer = {}

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

  const RENDER_DELAY = 250

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

  let currentRAF
  let tick = (t) => {
    t = currentTime + 1000 / fr
    currentTime = t

    app.tick(t)
    capturer.capture(app.canvas)

    if (currentTime <= 1000 * (secondsLong + captureTime) + 1000 / fr) {
      window.setTimeout(() => {
        if (!still) {
          currentRAF = window.requestAnimationFrame(tick)
        }
      }, RENDER_DELAY)
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
    if (!currentRAF) {
      currentRAF = window.requestAnimationFrame(tick)
    }
    run()
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
  gui.add(app, 'scale', -15, 15).step(0.0001).listen()
  gui.add(app, 'epsilon', 0.0000001, 0.05).step(0.0000001).listen()

  let angleCF = gui.addFolder('Angle Coefficients')
  angleCF.add(app, 'angle1C', -2, 2).step(0.00001).listen()
  angleCF.add(app, 'angle2C', -2, 2).step(0.001).listen()
  angleCF.add(app, 'angle3C', -2, 2).step(0.001).listen()

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
  cameraRotF.add(app, 'LOOKAT')
  cameraRotF.add(app.cameraAngles, '0', -Math.PI, Math.PI).step(0.001).listen()
  cameraRotF.add(app.cameraAngles, '1', -Math.PI, Math.PI).step(0.001).listen()
  cameraRotF.add(app.cameraAngles, '2', -Math.PI, Math.PI).step(0.001).listen()

  let colorsF = gui.addFolder('Colors')
  colorsF.addColor(app, 'colors1').listen()
  colorsF.addColor(app, 'colors2').listen()

  // ----- Setup -----
  window.addEventListener('resize', resize)

  // Run
  app.run()

  let currentRAF
  let tick = (t) => {
    currentRAF = window.requestAnimationFrame(tick)

    app.tick(t)
  }

  app.loaded.then(() => {
    if (!currentRAF) {
      currentRAF = window.requestAnimationFrame(tick)
    }
  })
}
