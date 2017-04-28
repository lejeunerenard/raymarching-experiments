const inquirer = require('inquirer')
const mkdirp = require('mkdirp')
const fs = require('fs')
const path = require('path')

const galleryDir = 'gallery'
const questions = [
  {
    type: 'input',
    name: 'name',
    default: () => {
      let today = new Date()
      return today.toISOString().replace(/T.*/, '')
    },
    message: 'Gallery Item Name'
  }
]
inquirer.prompt(questions).then((answers) => {
  let galleryItemDir = answers.name
  galleryItemDir = galleryItemDir.replace(/ /, '-').toLowerCase()
  galleryItemDir = path.join(galleryDir, galleryItemDir)

  mkdirp(galleryItemDir, (err) => {
    if (err && err.code !== 'EEXIST') {
      console.error(err)
      return
    }

    mkdirp(path.join(galleryItemDir, './node_modules/webvr-polyfill/build/'), (err) => {
      if (err && err.code !== 'EEXIST') {
        console.error(err)
        return
      }

      ['bundle.js', 'dat.gui.min.js', './node_modules/webvr-polyfill/build/webvr-polyfill.js', 'index.html'].forEach((file) => {
        fs.createReadStream(file)
          .pipe(fs.createWriteStream(path.join(galleryItemDir, file)))
      })
    })
  })
})
