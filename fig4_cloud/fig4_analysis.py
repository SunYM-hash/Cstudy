import os, glob, math
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from scipy.stats import kruskal

GROUPS=["CK","F-L20","C-L20","F-M20","C-M20","F-M40","C-M40"]
RAW_TO_GROUP={"CK":"CK","T2":"F-L20","T4":"C-L20","T6":"F-M20","T8":"C-M20","T5":"F-M40","T7":"C-M40"}
COLORS={"CK":"#1f77b4","F-L20":"#ff7f0e","C-L20":"#2ca02c","F-M20":"#d62728","C-M20":"#9467bd","F-M40":"#8c564b","C-M40":"#e377c2"}

def merge_parts(pattern):
    files=sorted(glob.glob(pattern))
    if not files:
        raise RuntimeError("No parts: "+pattern)
    out=None
    for fn in files:
        d=pd.read_csv(fn,sep="\t",index_col=0)
        out=d if out is None else out.add(d,fill_value=0)
    out=out.fillna(0)
    out=out.loc[out.sum(axis=1)>0]
    return out

def bh(p):
    p=np.asarray(p,float)
    n=len(p)
    order=np.argsort(p)
    q=np.empty(n,float)
    prev=1.0
    for rank_idx in range(n-1,-1,-1):
        i=order[rank_idx]
        rank=rank_idx+1
        val=min(prev,p[i]*n/rank)
        q[i]=val
        prev=val
    return np.clip(q,0,1)

def treatment_labels(cols):
    labs=[]
    for c in cols:
        raw=str(c).split("-")[0]
        labs.append(RAW_TO_GROUP[raw])
    return np.array(labs)

def analyze(counts, depth, min_prev=3, min_total=20):
    labs=treatment_labels(counts.columns)
    prev=(counts>0).sum(axis=1)
    tot=counts.sum(axis=1)
    x=counts.loc[(prev>=min_prev)&(tot>=min_total)].copy()

    # CLR on genus counts with small pseudocount; treatment tests on CLR values.
    logx=np.log(x+0.5)
    clr=logx.sub(logx.mean(axis=0),axis=1)

    # Relative abundance among classified genera, used only for intuitive heatmaps.
    rel=x.div(x.sum(axis=0),axis=1)*100
    gmean=pd.DataFrame(index=x.index,columns=GROUPS,dtype=float)
    for g in GROUPS:
        gmean[g]=rel.loc[:,labs==g].mean(axis=1)

    rows=[]
    for genus in x.index:
        vals=[clr.loc[genus,labs==g].values for g in GROUPS]
        try:
            H,p=kruskal(*vals)
        except Exception:
            H,p=np.nan,1.0
        eps=max(0.0,(H-len(GROUPS)+1)/(len(labs)-len(GROUPS))) if np.isfinite(H) else 0.0
        means=gmean.loc[genus]
        enriched=means.idxmax()
        max_ab=float(means.max())
        rows.append((genus,H,p,eps,enriched,max_ab))
    st=pd.DataFrame(rows,columns=["Genus","H","p","epsilon2","Enriched","MaxMeanRA"]).set_index("Genus")
    st["q"]=bh(st["p"].fillna(1).values)
    st["score"]=st["epsilon2"]*(-np.log10(st["p"].clip(lower=1e-12)))
    st["Depth"]=depth

    # candidate responsive genera: avoid extremely rare taxa.
    cand=st[(st["p"]<0.05)&(st["MaxMeanRA"]>=0.10)].copy()
    cand=cand.sort_values(["q","score","epsilon2"],ascending=[True,False,False])
    if len(cand)<8:
        # supplement only for visual preview; these are not called significant in caption.
        extra=st[(st["MaxMeanRA"]>=0.10)&(~st.index.isin(cand.index))].sort_values(["p","score"],ascending=[True,False])
        cand=pd.concat([cand,extra.head(8-len(cand))])
    sel=cand.head(12).index.tolist()

    hm=gmean.loc[sel].copy()
    z=hm.sub(hm.mean(axis=1),axis=0).div(hm.std(axis=1).replace(0,np.nan),axis=0).fillna(0)
    return st,gmean,z,sel

B20=merge_parts("fig4_cloud/parts/b20_part*.tsv")
B40=merge_parts("fig4_cloud/parts/b40_part*.tsv")
F20=merge_parts("fig4_cloud/parts/f20_part*.tsv")
F40=merge_parts("fig4_cloud/parts/f40_part*.tsv")

analyses={}
for key,mat,dep in [("b20",B20,"0–20 cm"),("b40",B40,"20–40 cm"),("f20",F20,"0–20 cm"),("f40",F40,"20–40 cm")]:
    analyses[key]=analyze(mat,dep)

out="fig4_cloud/output"
os.makedirs(out,exist_ok=True)
for key,(st,gm,z,sel) in analyses.items():
    st.to_csv(f"{out}/{key}_genus_stats.csv")
    gm.loc[sel].to_csv(f"{out}/{key}_selected_group_means.csv")

# Console summary
for key,(st,gm,z,sel) in analyses.items():
    print("\n",key,"tested",len(st),"p<0.05",int((st.p<0.05).sum()),"q<0.05",int((st.q<0.05).sum()))
    print("selected:",", ".join(sel))

# -------- figure --------
plt.rcParams.update({
    "font.family":"Liberation Serif",
    "font.serif":["Liberation Serif","Times New Roman","DejaVu Serif"],
    "mathtext.fontset":"stix",
    "font.size":10.0,
    "axes.labelsize":10.2,
    "xtick.labelsize":8.4,
    "ytick.labelsize":8.4,
    "axes.linewidth":0.9,
})
fig=plt.figure(figsize=(12.4,15.8),facecolor="white")
outer=fig.add_gridspec(3,2,left=.035,right=.985,top=.985,bottom=.025,wspace=.055,hspace=.095)

def panel_frame(cell,letter,title):
    fr=fig.add_subplot(cell)
    fr.set_xticks([]);fr.set_yticks([]);fr.set_facecolor("white")
    for sp in fr.spines.values():
        sp.set_visible(True);sp.set_color("black");sp.set_linewidth(1.0)
    fr.add_patch(Rectangle((0,.92),1,.08,transform=fr.transAxes,facecolor="#ededed",edgecolor="none"))
    fr.plot([0,1],[.92,.92],transform=fr.transAxes,color="#777",lw=.8,clip_on=False)
    fr.text(.028,.96,f"({letter})",transform=fr.transAxes,ha="left",va="center",fontsize=13,fontweight="bold")
    fr.text(.50,.96,title,transform=fr.transAxes,ha="center",va="center",fontsize=13.2,fontweight="bold")
    return fr

def add_rel(fr,rect):
    b=fr.get_position();l,bo,w,h=rect
    return fig.add_axes([b.x0+b.width*l,b.y0+b.height*bo,b.width*w,b.height*h])

def heat_panel(cell,letter,title,key):
    st,gm,z,sel=analyses[key]
    fr=panel_frame(cell,letter,title)
    ax=add_rel(fr,[.25,.16,.69,.70])
    im=ax.imshow(z.values,aspect="auto",cmap="RdBu_r",vmin=-2,vmax=2,interpolation="nearest")
    ax.set_xticks(np.arange(len(GROUPS)));ax.set_xticklabels(GROUPS,rotation=45,ha="right",fontsize=8.5)
    labels=[]
    for g in z.index:
        q=st.loc[g,"q"];p=st.loc[g,"p"]
        mark="**" if q<0.05 else ("*" if p<0.05 else "")
        labels.append(g+mark)
    ax.set_yticks(np.arange(len(labels)));ax.set_yticklabels(labels,fontstyle="italic",fontsize=8.2)
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(True);s.set_linewidth(.8);s.set_color("black")
    ax.set_xticks(np.arange(-.5,len(GROUPS),1),minor=True)
    ax.set_yticks(np.arange(-.5,len(labels),1),minor=True)
    ax.grid(which="minor",color="white",linewidth=.7)
    ax.tick_params(which="minor",bottom=False,left=False)
    cax=add_rel(fr,[.29,.065,.60,.025])
    cb=fig.colorbar(im,cax=cax,orientation="horizontal")
    cb.set_ticks([-2,0,2]);cb.ax.tick_params(labelsize=7.5,length=2)
    cb.set_label("Row Z-score of mean relative abundance",fontsize=8.3,labelpad=1)
    fr.text(.06,.095,"* p < 0.05; ** FDR q < 0.05",transform=fr.transAxes,fontsize=7.8,ha="left",va="center")

def summary_panel(cell,letter,title,keys):
    fr=panel_frame(cell,letter,title)
    dfs=[]
    for key in keys:
        st=analyses[key][0].copy()
        st=st[(st.p<0.05)&(st.MaxMeanRA>=0.10)].copy()
        st["Label"]=[f"{g} ({d})" for g,d in zip(st.index,st.Depth)]
        dfs.append(st)
    d=pd.concat(dfs) if dfs else pd.DataFrame()
    if len(d)==0:
        # fallback to strongest effects so panel is never empty
        d=pd.concat([analyses[k][0].head(0) for k in keys])
    d=d.sort_values(["q","score","epsilon2"],ascending=[True,False,False]).head(14)
    d=d.sort_values("epsilon2",ascending=True)
    ax=add_rel(fr,[.34,.14,.61,.72])
    y=np.arange(len(d))
    cols=[COLORS.get(g,"#777") for g in d["Enriched"]]
    ax.barh(y,d["epsilon2"].values,color=cols,height=.65)
    ax.set_yticks(y);ax.set_yticklabels(d["Label"].values,fontstyle="italic",fontsize=8.2)
    ax.set_xlabel(r"Treatment effect size ($\epsilon^2$)",fontweight="bold")
    ax.spines["top"].set_visible(False);ax.spines["right"].set_visible(False)
    ax.grid(axis="x",color=".88",lw=.6,zorder=0);ax.set_axisbelow(True)
    handles=[Line2D([0],[0],marker="s",linestyle="None",markersize=7,markerfacecolor=COLORS[g],markeredgecolor="none",label=g) for g in GROUPS]
    fr.legend(handles=handles,ncol=4,loc="lower center",bbox_to_anchor=(.5,.018),frameon=False,fontsize=7.7,handletextpad=.3,columnspacing=.7,title="Enriched treatment",title_fontsize=8.1)

heat_panel(outer[0,0],"A","Bacterial responsive genera — 0–20 cm","b20")
heat_panel(outer[0,1],"B","Bacterial responsive genera — 20–40 cm","b40")
heat_panel(outer[1,0],"C","Fungal responsive genera — 0–20 cm","f20")
heat_panel(outer[1,1],"D","Fungal responsive genera — 20–40 cm","f40")
summary_panel(outer[2,0],"E","Bacterial differential genera summary",["b20","b40"])
summary_panel(outer[2,1],"F","Fungal differential genera summary",["f20","f40"])

fig.savefig(f"{out}/Fig4_v1_preview.png",dpi=150,facecolor="white",bbox_inches="tight",pad_inches=.02)
fig.savefig(f"{out}/Fig4_v1_3x2.png",dpi=450,facecolor="white",bbox_inches="tight",pad_inches=.02)
fig.savefig(f"{out}/Fig4_v1_3x2.pdf",facecolor="white",bbox_inches="tight",pad_inches=.02)
fig.savefig(f"{out}/Fig4_v1_3x2_600dpi.tif",dpi=600,facecolor="white",bbox_inches="tight",pad_inches=.02,pil_kwargs={"compression":"tiff_lzw"})
plt.close(fig)
print("Fig4 complete")
