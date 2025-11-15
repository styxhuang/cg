async function j(url,opt){try{const r=await fetch(url,opt||{});const t=await r.text();try{return{ok:r.ok,status:r.status,data:JSON.parse(t)}}catch(e){return{ok:r.ok,status:r.status,data:{raw:t}}}}catch(e){return{ok:false,status:0,data:{error:String(e)}}}}
const rid=localStorage.getItem('rid')
const btnMap=document.getElementById('btnMap')
const btnRdf=document.getElementById('btnRdf')
const btnInvConfig=document.getElementById('btnInvConfig')
const btnInvRun=document.getElementById('btnInvRun')
const btnTable=document.getElementById('btnTable')
const btnShowCG=document.getElementById('btnShowCG')
const log=document.getElementById('log')
const vlog=document.getElementById('vlog')
function enable(el){if(el)el.disabled=false}
if(rid){[btnMap,btnRdf,btnInvConfig,btnInvRun,btnTable,btnShowCG].forEach(enable)}
btnMap.addEventListener('click',async()=>{const fd=new FormData();fd.append('rid',rid);const r=await j('/votca/map',{method:'POST',body:fd});log.textContent=JSON.stringify(r.data,null,2)})
btnRdf.addEventListener('click',async()=>{const fd=new FormData();fd.append('rid',rid);const r=await j('/rdf/run',{method:'POST',body:fd});log.textContent=JSON.stringify(r.data,null,2)})
btnInvConfig.addEventListener('click',async()=>{const fd=new FormData();fd.append('rid',rid);const r=await j('/inverse/config',{method:'POST',body:fd});log.textContent=JSON.stringify(r.data,null,2)})
btnInvRun.addEventListener('click',async()=>{const fd=new FormData();fd.append('rid',rid);const r=await j('/inverse/run',{method:'POST',body:fd});log.textContent=JSON.stringify(r.data,null,2)})
btnTable.addEventListener('click',async()=>{const fd=new FormData();fd.append('rid',rid);const r=await j('/ff/table',{method:'POST',body:fd});log.textContent=JSON.stringify(r.data,null,2)})
btnShowCG.addEventListener('click',async()=>{const r=await j('/file?path='+encodeURIComponent('backend/workspace/'+rid+'/cg_conf.gro'));const d=r.data;if(!d.ok){vlog.textContent=JSON.stringify(d);return}const conv=new FormData();conv.append('path','backend/workspace/'+rid+'/cg_conf.gro');const c=await j('/convert/pdb',{method:'POST',body:conv});if(!c.data.ok){vlog.textContent=JSON.stringify(c.data);return}const rd=await j('/file?path='+encodeURIComponent(c.data.out));const container=document.getElementById('viewer');container.innerHTML='';if(window.$3Dmol){const v=$3Dmol.createViewer('viewer',{backgroundColor:'white'});v.addModel(rd.data.content,'pdb');v.setStyle({}, {sphere:{radius:0.6}});v.zoomTo();v.render();vlog.textContent='已加载 CG'} })
