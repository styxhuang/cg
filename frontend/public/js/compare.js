async function j(url,opt){try{const r=await fetch(url,opt||{});const t=await r.text();try{return{ok:r.ok,status:r.status,data:JSON.parse(t)}}catch(e){return{ok:r.ok,status:r.status,data:{raw:t}}}}catch(e){return{ok:false,status:0,data:{error:String(e)}}}}
const rid=localStorage.getItem('rid')
const btnList=document.getElementById('btnList')
const btnCompare=document.getElementById('btnCompare')
const selComputed=document.getElementById('selComputed')
const target=document.getElementById('target')
const clog=document.getElementById('clog')
function enable(el){if(el)el.disabled=false}
if(rid){[btnList,btnCompare].forEach(enable)}
btnList.addEventListener('click',async()=>{const r=await j('/ff/list?rid='+encodeURIComponent(rid));selComputed.innerHTML='';(r.data.files||[]).forEach(f=>{const o=document.createElement('option');o.value=f;o.textContent=f;selComputed.appendChild(o)})})
btnCompare.addEventListener('click',async()=>{const comp=selComputed.value;const tgt=target.value;const fd1=new FormData();fd1.append('path',comp);const fd2=new FormData();fd2.append('path',tgt);const r1=await j('/ff/parse',{method:'POST',body:fd1});const r2=await j('/ff/parse',{method:'POST',body:fd2});clog.textContent='';const ds=[];if(r1.data.x&&r1.data.y)ds.push({x:r1.data.x,y:r1.data.y,type:'scatter',name:'Computed'});if(r2.data.x&&r2.data.y)ds.push({x:r2.data.x,y:r2.data.y,type:'scatter',name:'Target'});if(ds.length)Plotly.newPlot('plot',ds,{margin:{t:20},xaxis:{title:'r'},yaxis:{title:'g(r)'}})
})
