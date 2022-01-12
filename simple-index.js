import DefaultSceneRenderer from './default-scene-renderer'
import App from './app'

let app = new App({ sceneRenderer })
let sceneRenderer = new DefaultSceneRenderer(app.gl)
app.sceneRenderer = sceneRenderer.render.bind(sceneRenderer)
let resize = () => {
  let dim = app.getDimensions()
  sceneRenderer.resize(dim[0], dim[1])
}
window.addEventListener('resize', resize)
resize()

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
