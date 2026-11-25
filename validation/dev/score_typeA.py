# BPFEA_TESTS=<test data dir> python3 score_typeA.py <results dir>...: tall-model curves vs type A test means at fixed rotations, and RMS relative error per run tag
import csv,glob,sys,os,numpy as np
T=os.path.join(os.environ["BPFEA_TESTS"], "")   # folder with the test curves typeA_P<lb>_*.csv
th=[0.002,0.005,0.01,0.015,0.02,0.03,0.04,0.06,0.08]
def test(N):
    out=[]
    for f in sorted(glob.glob(T+f"typeA_P{N*1000}_*")):
        a=np.genfromtxt(f,delimiter=',',skip_header=1)[:,:2]; a=a[~np.isnan(a).any(1)]
        i=np.argmax(a[:,0]); a=a[:i+1]
        out.append([np.interp(t,a[:,0],a[:,1]) if t<=a[-1,0] else np.nan for t in th])
    return np.nanmean(out,0)
def model(f,col='M_kipin'):
    r=list(csv.DictReader(open(f))); x=np.array([float(q['rotation_rad']) for q in r]); y=np.array([float(q[col]) for q in r])
    return np.array([np.interp(t,x,y) if t<=x.max() else np.nan for t in th])
P=lambda l,v: print("%-30s"%l,*["%6.1f"%q for q in v])
err={}
for N in (5,10,15):
    print("%-30s"%f"N={N}  theta",*["%6.3f"%t for t in th]); t=test(N); P("test mean",t)
    for d in sys.argv[1:]:
        for f in sorted(glob.glob(f"{d}/A_N{N}_*.csv")):
            tag=os.path.basename(f)[len(f"A_N{N}_"):-4]; m=model(f); P(tag,m); P("    anchor kip",model(f,'bolt1_kip'))
            err.setdefault(tag,{})[N]=np.sqrt(np.nanmean((m/t-1)**2)) if not np.isnan(m).all() else np.nan
print("\nRMS relative error    N=5   N=10   N=15   all")
for tag,e in sorted(err.items(), key=lambda kv: np.nanmean(list(kv[1].values()))):
    v=[e.get(N,np.nan) for N in (5,10,15)]; print("%-28s"%tag,*["%5.1f%%"%(100*x) for x in v],"%5.1f%%"%(100*np.sqrt(np.nanmean(np.square(v)))))
