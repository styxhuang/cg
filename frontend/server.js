const express = require('express')
const path = require('path')
const layouts = require('express-ejs-layouts')
const { createProxyMiddleware } = require('http-proxy-middleware')
const app = express()
const PORT = process.env.PORT || 5173
const BACKEND = process.env.BACKEND || 'http://localhost:5050'

app.set('view engine','ejs')
app.set('views', path.join(__dirname, 'views'))
app.use(layouts)
app.set('layout', 'layout')
app.use(express.static(path.join(__dirname, 'public')))
app.use((req,res,next)=>{ try { console.log(`[frontend] ${req.method} ${req.originalUrl}`) } catch(e){} next() })
const API_PATHS = ['/upload','/env','/status','/votca','/workflow','/file','/ff','/inverse','/rdf','/cg','/convert','/routes']
const apiProxyOpts = {
  target: BACKEND,
  changeOrigin: true,
  logLevel: 'debug',
  onError(err, req, res){ try { console.error(`[proxy api] error`, err && err.message) } catch(e){} },
  onProxyReq(proxyReq, req, res){ try { console.log(`[proxy api] req ${req.method} ${req.originalUrl} -> ${BACKEND}`) } catch(e){} },
  onProxyRes(proxyRes, req, res){ try { console.log(`[proxy api] res ${proxyRes.statusCode} for ${req.originalUrl}`) } catch(e){} }
}
console.log(`[proxy] mount contexts -> ${BACKEND}`)
const apiProxy = createProxyMiddleware((pathname, req) => API_PATHS.some(p => pathname.startsWith(p)), apiProxyOpts)
app.use(apiProxy)

app.get('/', (req,res)=>{ res.render('index', { backend: BACKEND, title: '上传' }) })
app.get('/workflow', (req,res)=>{ res.render('workflow', { backend: BACKEND, title: '工作流' }) })
app.get('/compare', (req,res)=>{ res.render('compare', { backend: BACKEND, title: '对比' }) })

app.listen(PORT, ()=>{ process.stdout.write(`http://localhost:${PORT}/`) })
