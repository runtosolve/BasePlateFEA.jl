# Build the validation report page: moment–rotation curves of the Barcelona tests, BasePlateFEA
# (elastic and elastic–plastic) and S. Ádány's MATLAB code, one panel per plate type and axial load.
#   julia --project=. validation/make_report.jl <validation dir>
using DelimitedFiles, Printf
import JSON3
include(joinpath(@__DIR__, "barcelona_cases.jl"))

interp(x, xs, ys) = (i = clamp(searchsortedlast(xs, x), 1, length(xs) - 1); ys[i] + (ys[i+1] - ys[i]) * (x - xs[i]) / (xs[i+1] - xs[i]))
r5(v) = round.(v; sigdigits = 5)

function model_curve(path)
    isfile(path) || return nothing
    d = readdlm(path, ',', Float64; skipstart = 1)
    d = d[findlast(==(1), d[:, 1]):end, :]
    return (θ = d[:, 4], M = d[:, 3], peeq = d[:, 9], bolt = max.(d[:, 6], d[:, 7]), contact = d[:, 8])
end

function main(vdir)
    resdir = joinpath(vdir, "results"); datadir = joinpath(vdir, "barcelona_test_data"); octdir = joinpath(vdir, "octave_reference")
    panels = []; rows = []
    for typ in (:A, :B, :C), N in BARCELONA[typ].loads
        tag = "$(typ)_N$(loadtag(N))"
        tests = []
        for f in sort(readdir(datadir))
            startswith(f, "type$(typ)_P$(round(Int, 1000N))_") && endswith(f, ".csv") || continue
            d = readdlm(joinpath(datadir, f), ',', Float64; skipstart = 1)
            n = argmax(d[:, 1])
            push!(tests, (θ = d[1:n, 1], M = d[1:n, 2]))
        end
        pl = model_curve(joinpath(resdir, tag * "_plastic.csv")); el = model_curve(joinpath(resdir, tag * "_elastic.csv"))
        oct = nothing
        op = joinpath(octdir, "$(typ)_N$(loadtag(N))_m_r.txt")
        if isfile(op)
            o = readdlm(op, ','); o = o[findlast(==(0), o[:, 1]):end, :]
            oct = (θ = o[:, 2], M = o[:, 1])
        end
        tall = nothing                     # type A with the upright modelled to the load point (run_typeA_tall.jl)
        tp = joinpath(vdir, "results_tall", tag * "_tall_imp03.csv")
        if isfile(tp)
            d = readdlm(tp, ',', Float64; skipstart = 1); d = d[findlast(==(1), d[:, 1]):end, :]
            tall = (θ = d[:, 5], M = d[:, 4])
        end
        series = Any[]
        for t in tests
            push!(series, Dict(:kind => "test", :x => r5(t.θ), :y => r5(t.M)))
        end
        oct === nothing || push!(series, Dict(:kind => "matlab", :x => r5(oct.θ), :y => r5(oct.M)))
        el === nothing || push!(series, Dict(:kind => "elastic", :x => r5(el.θ), :y => r5(el.M)))
        pl === nothing || push!(series, Dict(:kind => "plastic", :x => r5(pl.θ), :y => r5(pl.M)))
        tall === nothing || push!(series, Dict(:kind => "tall", :x => r5(tall.θ), :y => r5(tall.M)))
        push!(panels, Dict(:type => String(typ), :N => N, :tag => tag, :series => series, :xmax => BARCELONA[typ].θmax,
                           :viewer => isfile(joinpath(vdir, "viewers_makie", tag * "_plastic.html")) ? "viewers_makie/$(tag)_plastic.html" : nothing,
                           :viewer_tall => isfile(joinpath(vdir, "viewers_makie", tag * "_tall_imp03.html")) ? "viewers_makie/$(tag)_tall_imp03.html" : nothing))
        # moments at fixed rotations: model / mean of the tests that reach the rotation
        if pl !== nothing && !isempty(tests)
            for θ in (0.002, 0.005, 0.01, 0.02, 0.04)
                θ <= pl.θ[end] || continue
                tm = [interp(θ, monotone(t.θ, t.M)...) for t in tests if maximum(t.θ) >= θ]
                isempty(tm) && continue
                Mt = sum(tm) / length(tm)
                push!(rows, Dict(:type => String(typ), :N => N, :theta => θ, :ntests => length(tm), :Mtest => round(Mt, digits = 2),
                                 :Mplastic => round(interp(θ, pl.θ, pl.M), digits = 2),
                                 :Mtall => tall !== nothing && θ <= tall.θ[end] ? round(interp(θ, tall.θ, tall.M), digits = 2) : nothing,
                                 :Melastic => el !== nothing && θ <= el.θ[end] ? round(interp(θ, el.θ, el.M), digits = 2) : nothing))
            end
        end
    end
    html = replace(TEMPLATE, "/*DATA*/null" => replace(JSON3.write(Dict(:panels => panels, :rows => rows)), "</" => "<\\/"))
    p = joinpath(vdir, "validation_report.html")
    write(p, html)
    println(p)
    # console summary
    for typ in ("A", "B", "C"), θ in (0.005, 0.02, 0.04)
        rs = [r for r in rows if r[:type] == typ && r[:theta] == θ]
        isempty(rs) && continue
        q = [r[:Mplastic] / r[:Mtest] for r in rs]
        @printf("type %s  θ = %.3f  M_model/M_test: mean %.2f  range %.2f – %.2f  (%d load levels)\n", typ, θ, sum(q) / length(q), minimum(q), maximum(q), length(q))
    end
end

"running maximum of the rotation so that the test curve can be interpolated"
function monotone(θ, M)
    x = Float64[]; y = Float64[]
    for i in eachindex(θ)
        (isempty(x) || θ[i] > x[end]) && (push!(x, θ[i]); push!(y, M[i]))
    end
    return x, y
end

const TEMPLATE = raw"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>BasePlateFEA Validation</title>
<style>
:root { --page:#f9f9f7; --surface:#fcfcfb; --ink:#0b0b0b; --ink2:#52514e; --muted:#898781; --grid:#e1e0d9; --axis:#c3c2b7; --s1:#2a78d6; --s2:#eb6834; --s3:#1baf7a; --s4:#8a4fd0; }
@media (prefers-color-scheme: dark) { :root:not([data-theme="light"]) { --page:#0d0d0d; --surface:#1a1a19; --ink:#f0efec; --ink2:#c3c2b7; --grid:#2c2c2a; --axis:#383835; --s1:#3987e5; --s2:#d95926; --s3:#199e70; --s4:#a678e8; } }
:root[data-theme="dark"] { --page:#0d0d0d; --surface:#1a1a19; --ink:#f0efec; --ink2:#c3c2b7; --grid:#2c2c2a; --axis:#383835; --s1:#3987e5; --s2:#d95926; --s3:#199e70; --s4:#a678e8; }
body { margin:0; background:var(--page); color:var(--ink); font:14px/1.5 system-ui,-apple-system,"Segoe UI",sans-serif; }
main { max-width:1180px; margin:0 auto; padding:24px 16px 48px; }
h1 { font-size:22px; margin:0 0 4px; } h2 { font-size:16px; margin:28px 0 6px; } p { color:var(--ink2); max-width:78ch; margin:6px 0; }
.key { display:flex; gap:18px; flex-wrap:wrap; color:var(--ink2); margin:10px 0; }
.key i { display:inline-block; width:22px; border-top:2px solid; vertical-align:middle; margin-right:6px; }
.grid { display:grid; grid-template-columns:repeat(auto-fill,minmax(340px,1fr)); gap:12px; }
.card { background:var(--surface); border:1px solid var(--grid); border-radius:8px; padding:10px 12px; }
.card h3 { font-size:13px; margin:0 0 2px; display:flex; justify-content:space-between; font-weight:600; }
.card a { color:var(--s1); font-weight:400; }
svg { width:100%; height:auto; display:block; } svg text { fill:var(--muted); font-size:11px; }
.tw { overflow-x:auto; } table { border-collapse:collapse; font-variant-numeric:tabular-nums; background:var(--surface); }
td,th { text-align:right; padding:3px 10px; border-bottom:1px solid var(--grid); } th { color:var(--ink2); font-weight:500; }
#tip { position:fixed; pointer-events:none; background:var(--surface); border:1px solid var(--axis); border-radius:4px; padding:4px 8px; display:none; font-size:12px; }
</style></head><body><main>
<h1>BasePlateFEA validation: Barcelona base plate tests</h1>
<p>Moment–rotation response of base plate types A, B and C under constant axial compression. Moment in kip-in about the upright's axis of symmetry, rotation in rad between the two flange measurement points.</p>
<div class="key">
<span><i style="border-color:var(--s1)"></i>BasePlateFEA, elastic–plastic steel</span>
<span><i style="border-color:var(--s4);border-top-style:dashed"></i>BasePlateFEA, type A with 24 in upright, 0.3t imperfection, Yun–Gardner steel</span>
<span><i style="border-color:var(--s2);border-top-style:dashed"></i>BasePlateFEA, elastic steel</span>
<span><i style="border-color:var(--s3);border-top-style:dotted"></i>Ádány MATLAB code (elastic, nominal thickness)</span>
<span><i style="border-color:var(--muted);border-top-width:1px"></i>Tests (each repeat)</span></div>
<div id="sections"></div>
<h2>Model moment relative to the tests at fixed rotations</h2>
<p>Test moment is the mean of the repeats that reach the rotation. Test rotations are used as recorded (not zeroed).</p>
<div class="tw"><table id="tbl"></table></div>
</main><div id="tip"></div>
<script>
const D = /*DATA*/null;
const ST = { test:['--muted',1,''], matlab:['--s3',1.8,'2 3'], elastic:['--s2',1.8,'6 3'], plastic:['--s1',2.2,''], tall:['--s4',2.2,'9 3'] };
const NAME = { test:'test', matlab:'MATLAB', elastic:'elastic', plastic:'elastic–plastic', tall:'tall upright' };
function nice(v){const p=Math.pow(10,Math.floor(Math.log10(v)));const m=v/p;return (m<=1?1:m<=2?2:m<=5?5:10)*p;}
function chart(el,p){
  const W=360,H=230,m={l:38,r:14,t:8,b:26};
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
    const seen={}; let t=`θ = ${th.toFixed(4)} rad`;
    p.series.forEach(s=>{ if(s.kind==='test'||seen[s.kind])return; seen[s.kind]=1; let k=s.x.findIndex(x=>x>=th); if(k<1)return;
      const y=s.y[k-1]+(s.y[k]-s.y[k-1])*(th-s.x[k-1])/(s.x[k]-s.x[k-1]); t+=`<br>${NAME[s.kind]}: ${y.toFixed(1)} kip-in`; });
    tip.innerHTML=t; tip.style.display='block'; tip.style.left=e.clientX+12+'px'; tip.style.top=e.clientY+12+'px';};
  svg.onmouseleave=()=>tip.style.display='none';
}
const sec=document.getElementById('sections');
['A','B','C'].forEach(t=>{ const ps=D.panels.filter(p=>p.type===t); if(!ps.length)return;
  sec.insertAdjacentHTML('beforeend',`<h2>Type ${t}</h2><div class="grid" id="g${t}"></div>`);
  ps.forEach(p=>{ const c=document.createElement('div'); c.className='card';
    c.innerHTML=`<h3><span>Type ${t}, N = ${p.N} kip</span><span>${p.viewer?`<a href="${p.viewer}">3D viewer</a>`:''}${p.viewer_tall?` <a href="${p.viewer_tall}">3D, tall upright</a>`:''}</span></h3><div></div>`;
    document.getElementById('g'+t).appendChild(c); chart(c.lastChild,p); }); });
let h='<tr><th>Type</th><th>N (kip)</th><th>θ (rad)</th><th>Tests</th><th>M test</th><th>M elastic–plastic</th><th>ratio</th><th>M elastic</th><th>ratio</th><th>M tall upright</th><th>ratio</th></tr>';
D.rows.forEach(r=>h+=`<tr><td>${r.type}</td><td>${r.N}</td><td>${r.theta}</td><td>${r.ntests}</td><td>${r.Mtest.toFixed(2)}</td><td>${r.Mplastic.toFixed(2)}</td><td>${(r.Mplastic/r.Mtest).toFixed(2)}</td><td>${r.Melastic==null?'–':r.Melastic.toFixed(2)}</td><td>${r.Melastic==null?'–':(r.Melastic/r.Mtest).toFixed(2)}</td><td>${r.Mtall==null?'–':r.Mtall.toFixed(2)}</td><td>${r.Mtall==null?'–':(r.Mtall/r.Mtest).toFixed(2)}</td></tr>`);
document.getElementById('tbl').innerHTML=h;
</script></body></html>
"""

main(ARGS[1])
