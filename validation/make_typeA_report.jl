# Type A study page: tests against the 3 in stub model and the tall-upright models
# (run_typeA_tall.jl): imperfection amplitude and anchor stiffness variants.
#   julia --project=. validation/make_typeA_report.jl <validation dir>
using DelimitedFiles, Printf
import JSON3
include(joinpath(@__DIR__, "barcelona_cases.jl"))

interp(x, xs, ys) = (i = clamp(searchsortedlast(xs, x), 1, length(xs) - 1); ys[i] + (ys[i+1] - ys[i]) * (x - xs[i]) / (xs[i+1] - xs[i]))
r5(v) = round.(v; sigdigits = 5)
function monotone(θ, M)
    x = Float64[]; y = Float64[]
    for i in eachindex(θ)
        (isempty(x) || θ[i] > x[end]) && (push!(x, θ[i]); push!(y, M[i]))
    end
    return x, y
end
function curve(path, cθ, cM)
    isfile(path) || return nothing
    d = readdlm(path, ',', Float64; skipstart = 1)
    d = d[findlast(==(1), d[:, 1]):end, :]
    return (θ = d[:, cθ], M = d[:, cM])
end

const VARIANTS = [("stub", "3 in stub (previous model)"), ("tall_imp0", "tall upright, no imperfection"), ("tall_imp03", "tall, imperfection 0.3 t"),
                  ("tall_imp10", "tall, imperfection 1.0 t"), ("tall_imp03_bolt685", "tall, 0.3 t, anchors 685 kip/in")]

function main(vdir)
    datadir = joinpath(vdir, "barcelona_test_data")
    panels = []; rows = []
    for N in BARCELONA[:A].loads
        tag = "A_N$(loadtag(N))"
        tests = []
        for f in sort(readdir(datadir))
            startswith(f, "typeA_P$(round(Int, 1000N))_") && endswith(f, ".csv") || continue
            d = readdlm(joinpath(datadir, f), ',', Float64; skipstart = 1)
            n = argmax(d[:, 1])
            push!(tests, (θ = d[1:n, 1], M = d[1:n, 2]))
        end
        series = Any[Dict(:kind => "test", :x => r5(t.θ), :y => r5(t.M)) for t in tests]
        for (v, _) in VARIANTS
            c = v == "stub" ? curve(joinpath(vdir, "results", tag * "_plastic.csv"), 4, 3) : curve(joinpath(vdir, "results_tall", "$(tag)_$(v).csv"), 5, 4)
            c === nothing && continue
            push!(series, Dict(:kind => v, :x => r5(c.θ), :y => r5(c.M)))
            for θ in (0.005, 0.01, 0.02, 0.03, 0.04, 0.06)
                θ <= c.θ[end] || continue
                tm = [interp(θ, monotone(t.θ, t.M)...) for t in tests if maximum(t.θ) >= θ]
                isempty(tm) && continue
                push!(rows, Dict(:N => N, :theta => θ, :variant => v, :Mtest => sum(tm) / length(tm), :M => interp(θ, c.θ, c.M)))
            end
        end
        push!(panels, Dict(:N => N, :tag => tag, :series => series, :xmax => 0.08))
    end
    data = Dict(:panels => panels, :rows => rows, :variants => [Dict(:key => k, :label => l) for (k, l) in VARIANTS])
    p = joinpath(vdir, "typeA_study.html")
    write(p, replace(TEMPLATE, "/*DATA*/null" => replace(JSON3.write(data), "</" => "<\\/")))
    println(p)
    for (v, l) in VARIANTS, θ in (0.005, 0.02, 0.04)
        q = [r[:M] / r[:Mtest] for r in rows if r[:variant] == v && r[:theta] == θ]
        isempty(q) || @printf("%-34s θ = %.3f  M_model/M_test: %s\n", l, θ, join([@sprintf("%.2f", x) for x in q], "  "))
    end
end

const TEMPLATE = raw"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>Type A Study</title>
<style>
:root { --page:#f9f9f7; --surface:#fcfcfb; --ink:#0b0b0b; --ink2:#52514e; --muted:#898781; --grid:#e1e0d9; --axis:#c3c2b7; --b1:#86b6ef; --b2:#2a78d6; --b3:#0d366b; --s2:#eb6834; --s3:#1baf7a; }
@media (prefers-color-scheme: dark) { :root:not([data-theme="light"]) { --page:#0d0d0d; --surface:#1a1a19; --ink:#f0efec; --ink2:#c3c2b7; --grid:#2c2c2a; --axis:#383835; --b1:#1c5cab; --b2:#3987e5; --b3:#b7d3f6; --s2:#d95926; --s3:#199e70; } }
:root[data-theme="dark"] { --page:#0d0d0d; --surface:#1a1a19; --ink:#f0efec; --ink2:#c3c2b7; --grid:#2c2c2a; --axis:#383835; --b1:#1c5cab; --b2:#3987e5; --b3:#b7d3f6; --s2:#d95926; --s3:#199e70; }
body { margin:0; background:var(--page); color:var(--ink); font:14px/1.5 system-ui,-apple-system,"Segoe UI",sans-serif; }
main { max-width:1180px; margin:0 auto; padding:24px 16px 48px; }
h1 { font-size:22px; margin:0 0 4px; } h2 { font-size:16px; margin:28px 0 6px; } p, li { color:var(--ink2); max-width:78ch; } p { margin:6px 0; }
.key { display:flex; gap:6px 18px; flex-wrap:wrap; color:var(--ink2); margin:10px 0; }
.key i { display:inline-block; width:22px; border-top:2px solid; vertical-align:middle; margin-right:6px; }
.grid { display:grid; grid-template-columns:repeat(auto-fill,minmax(320px,1fr)); gap:12px; }
.card { background:var(--surface); border:1px solid var(--grid); border-radius:8px; padding:10px 12px; }
.card h3 { font-size:13px; margin:0 0 2px; font-weight:600; }
svg { width:100%; height:auto; display:block; } svg text { fill:var(--muted); font-size:11px; }
.tw { overflow-x:auto; } table { border-collapse:collapse; font-variant-numeric:tabular-nums; background:var(--surface); }
td,th { text-align:right; padding:3px 10px; border-bottom:1px solid var(--grid); white-space:nowrap; } th { color:var(--ink2); font-weight:500; } td:first-child,th:first-child { text-align:left; }
#tip { position:fixed; pointer-events:none; background:var(--surface); border:1px solid var(--axis); border-radius:4px; padding:4px 8px; display:none; font-size:12px; }
</style></head><body><main>
<h1>Type A: tall upright, 85 mm rotation, imperfection</h1>
<p>Moment–rotation response of Barcelona base plate type A. Moment in kip-in at the base plate; for the tall-upright models it includes the axial force times the lateral displacement of the loaded end, as in the test definition. Rotation in rad.</p>
<div class="key" id="key"></div>
<div class="grid" id="grid"></div>
<h2>Model moment divided by mean test moment</h2>
<p>Test moment is the mean of the repeats that reach the rotation; test rotations as recorded. A dash means the model run did not reach that rotation.</p>
<div class="tw"><table id="tbl"></table></div>
</main><div id="tip"></div>
<script>
const D = /*DATA*/null;
const ST = { test:['--muted',1,''], stub:['--s2',1.8,'6 3'], tall_imp0:['--b1',2,''], tall_imp03:['--b2',2.4,''], tall_imp10:['--b3',2,''], tall_imp03_bolt685:['--s3',2,'2 3'] };
const NAME = Object.fromEntries(D.variants.map(v=>[v.key,v.label])); NAME.test='Tests (each repeat)';
document.getElementById('key').innerHTML=['test',...D.variants.map(v=>v.key)].map(k=>`<span><i style="border-color:var(${ST[k][0]});border-top-width:${Math.min(ST[k][1],2)}px;border-top-style:${ST[k][2]==='6 3'?'dashed':ST[k][2]?'dotted':'solid'}"></i>${NAME[k]}</span>`).join('');
function nice(v){const p=Math.pow(10,Math.floor(Math.log10(v)));const m=v/p;return (m<=1?1:m<=2?2:m<=5?5:10)*p;}
function chart(el,p){
  const W=360,H=250,m={l:38,r:14,t:8,b:26};
  const xm=p.xmax, ym=nice(Math.max(...p.series.flatMap(s=>s.y.filter((y,i)=>s.x[i]<=xm)))*1.02);
  const X=v=>m.l+(W-m.l-m.r)*v/xm, Y=v=>H-m.b-(H-m.t-m.b)*v/ym; let h='';
  for(let i=0;i<=4;i++){const gx=xm*i/4,gy=ym*i/4;
    h+=`<line x1="${X(0)}" x2="${X(xm)}" y1="${Y(gy)}" y2="${Y(gy)}" stroke="var(${i?'--grid':'--axis'})"/>`;
    h+=`<text x="${X(0)-5}" y="${Y(gy)+4}" text-anchor="end">${+gy.toPrecision(3)}</text><text x="${X(gx)}" y="${H-m.b+14}" text-anchor="middle">${+gx.toPrecision(3)}</text>`;}
  h+=`<clipPath id="c${p.tag}"><rect x="${X(0)}" y="${m.t}" width="${X(xm)-X(0)}" height="${Y(0)-m.t}"/></clipPath>`;
  p.series.forEach(s=>{const st=ST[s.kind];
    h+=`<polyline clip-path="url(#c${p.tag})" fill="none" stroke="var(${st[0]})" stroke-width="${st[1]}" stroke-dasharray="${st[2]}" stroke-linejoin="round" points="${s.x.map((x,i)=>X(x).toFixed(1)+','+Y(s.y[i]).toFixed(1)).join(' ')}"/>`;});
  el.innerHTML=`<svg viewBox="0 0 ${W} ${H}">${h}</svg>`;
  const svg=el.firstChild, tip=document.getElementById('tip');
  svg.onmousemove=e=>{const r=svg.getBoundingClientRect(),th=((e.clientX-r.left)/r.width*W-m.l)/(W-m.l-m.r)*xm; if(th<0||th>xm){tip.style.display='none';return;}
    let t=`θ = ${th.toFixed(4)} rad`;
    p.series.forEach(s=>{ if(s.kind==='test')return; let k=s.x.findIndex(x=>x>=th); if(k<1)return;
      const y=s.y[k-1]+(s.y[k]-s.y[k-1])*(th-s.x[k-1])/(s.x[k]-s.x[k-1]); t+=`<br>${NAME[s.kind]}: ${y.toFixed(1)}`; });
    tip.innerHTML=t; tip.style.display='block'; tip.style.left=Math.min(e.clientX+12,innerWidth-250)+'px'; tip.style.top=e.clientY+12+'px';};
  svg.onmouseleave=()=>tip.style.display='none';
}
D.panels.forEach(p=>{ const c=document.createElement('div'); c.className='card';
  c.innerHTML=`<h3>Type A, N = ${p.N} kip</h3><div></div>`; document.getElementById('grid').appendChild(c); chart(c.lastChild,p); });
const Ns=[...new Set(D.rows.map(r=>r.N))], ths=[...new Set(D.rows.map(r=>r.theta))].sort((a,b)=>a-b);
let h='<tr><th>Model</th><th>N (kip)</th>'+ths.map(t=>`<th>θ = ${t}</th>`).join('')+'</tr>';
D.variants.forEach(v=>Ns.forEach(N=>{ const rs=D.rows.filter(r=>r.variant===v.key&&r.N===N); if(!rs.length)return;
  h+=`<tr><td>${v.label}</td><td>${N}</td>`+ths.map(t=>{const r=rs.find(r=>r.theta===t);return `<td>${r?(r.M/r.Mtest).toFixed(2):'–'}</td>`}).join('')+'</tr>';}));
document.getElementById('tbl').innerHTML=h;
</script></body></html>
"""
main(ARGS[1])
