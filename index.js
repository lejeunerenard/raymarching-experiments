import DefaultSceneRenderer from './default-scene-renderer'
import CCaptureCapturer from './capturer/ccapture'
import App from './app'

import { name } from './info.json'

const fr = 60
const captureTime = 0 * 5
const secondsLong = 15
const capturing = false

let app = new App()
window.app = app

app.width = window.devicePixelRatio * 1080
app.height = window.devicePixelRatio * 1920

let sceneRenderer = new DefaultSceneRenderer(app.gl)
app.sceneRenderer = sceneRenderer.render.bind(sceneRenderer)

// Setup resizing
let resize = () => {
  let dim = app.getDimensions()
  sceneRenderer.resize(dim[0], dim[1])
}
resize()

app.totalTime = secondsLong
// window.time = 0.3
const still = false

if (capturing) {
  let massagedName = name.replace(/ /g, '-')
  massagedName = massagedName.replace(/'/g, '')
  massagedName = massagedName.toLowerCase()

  let renderCount = 3

  let renderings = {}
  try {
    renderings = JSON.parse(window.localStorage.renderings)
  } catch (e) {}
  if (!renderings) renderings = {}

  if (massagedName in renderings) {
    renderCount = renderings[massagedName] + 1
  }

  renderings[massagedName] = renderCount

  let filename = massagedName + '-render' + renderCount
  console.log('filename', filename)
  let capturer = new CCaptureCapturer({
    format: 'jpg',
    framerate: fr,
    filename,
    autoSaveTime: 5,
    quality: 98,
    startTime: captureTime,
    timeLength: secondsLong
  })

  let currentTime = captureTime * 1000
  window.capturer = capturer

  let run = () => {
    capturer.start()
    app.run()
  }

  let currentRAF
  let tick = async (t) => {
    t = currentTime + 1000 / fr
    currentTime = t

    app.tick(t)
    await capturer.capture(app.canvas)

    if (currentTime <= 1000 * (secondsLong + captureTime) + 1000 / fr) {
      if (!still) {
        currentRAF = window.requestAnimationFrame(tick)
      }
    } else {
      capturer.save().then(() => {
        window.localStorage.renderings = JSON.stringify(renderings)
        window.setTimeout(() => {
          console.log('done sending message')
          var xmlHTTP = new XMLHttpRequest()
          xmlHTTP.open('GET', 'http://' + window.location.hostname + ':7321/', false)
          xmlHTTP.send(null)
          console.log('response text', xmlHTTP.responseText)
        }, 250)
      })
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
  angleCF.add(app, 'angle1C', -5, 5).step(0.0001).listen()
  angleCF.add(app, 'angle2C', -5, 5).step(0.0001).listen()
  angleCF.add(app, 'angle3C', -5, 5).step(0.001).listen()

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
