import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

D = {"g":["CK","F-L20","C-L20","F-M20","C-M20","F-M40","C-M40"],"bt":["Proteobacteria","Acidobacteriota","Actinobacteriota","Chloroflexi","Gemmatimonadota","Myxococcota","Bacteroidota","Other taxa"],"ft":["Ascomycota","Basidiomycota","Mortierellomycota","Glomeromycota","Chytridiomycota","Rozellomycota","Other taxa"],"b20":[[26.91,26.14,20.62,32.05,26.6,27.7,29.07],[24.31,23.67,27.35,23.05,21.58,20.75,23.7],[16.07,13.17,11.3,10.28,9.94,12.81,11.64],[7.97,7.52,7.58,5.65,6.27,5.33,7.17],[4.4,5.93,8.13,4.07,4.1,5.28,5.4],[3.11,4.32,3.89,4.79,3.81,5.13,4.13],[2.51,2.2,1.63,3.7,4.49,2.99,3.97],[14.71,17.06,19.5,16.4,23.21,19.99,14.91]],"b40":[[24.51,25.62,32.02,30.13,26.84,37.65,29.6],[24.3,23.52,13.61,22.9,23.45,17.18,17.55],[14.32,10.7,13.17,9.87,9.98,9.99,13.68],[8.11,7.5,10.77,6.22,9.15,4.26,8.2],[6.26,6.05,5.88,5.68,5.95,5.65,6.14],[3.41,3.81,3.42,3.76,3.72,3.71,4.37],[1.79,1.78,5.91,2.72,3.53,3.69,4.94],[17.31,21.04,15.24,18.71,17.38,17.88,15.51]],"f20":[[21.08,27.24,33.91,24.6,36.68,44.79,40.32],[68.29,46.73,27.54,48.03,49.32,21.41,32.65],[5.85,14.6,23.66,14.88,9.04,22.94,10.33],[1.17,2.45,2.44,1.53,1.05,1.46,1.58],[0.11,1.16,2.25,2.64,0.99,1.22,7.52],[0.5,0.42,1.23,0.32,0.17,0.69,0.82],[3.01,7.4,8.97,8.02,2.74,7.49,6.78]],"f40":[[57.63,42.94,57.83,47.62,51.46,36.24,34.46],[25.74,19.82,18.18,30.67,21.31,44.44,37.55],[5.32,24.93,11.64,10.54,18.78,13.63,9.17],[0.85,2.32,1.38,1.3,2.78,0.94,2.89],[0.68,0.73,0.73,2.04,0.62,0.57,1.95],[0.35,0.46,1.26,1.07,0.69,0.85,0.65],[9.43,8.8,8.98,6.76,4.35,3.32,13.33]],"bs":[[[7.121,0.025],[7.182,0.03],[7.099,0.089],[7.2,0.036],[6.922,0.023],[7.232,0.01],[7.234,0.018]],[[7.157,0.061],[7.066,0.035],[7.025,0.013],[7.158,0.04],[7.239,0.079],[6.973,0.065],[7.197,0.028]]],"fs":[[[3.515,0.143],[4.076,0.055],[4.029,0.13],[4.186,0.137],[3.966,0.015],[3.937,0.061],[4.734,0.025]],[[4.923,0.027],[4.449,0.034],[4.812,0.087],[5.437,0.003],[5.021,0.022],[3.363,0.071],[4.226,0.012]]],"bp20":[[0,-0.093,0.109],[0,-0.09,0.134],[0,-0.093,0.141],[1,-0.111,0.013],[1,-0.172,-0.03],[1,-0.084,0.048],[2,-0.239,-0.249],[2,-0.211,-0.219],[2,-0.204,-0.208],[5,-0.076,-0.053],[5,-0.069,-0.04],[5,-0.047,-0.032],[3,-0.006,0.181],[3,-0.019,0.16],[3,-0.028,0.16],[6,0.217,0.139],[6,0.175,0.156],[6,0.113,0.115],[4,0.334,-0.165],[4,0.331,-0.159],[4,0.373,-0.202]],"bp40":[[0,0.205,-0.069],[0,0.203,-0.072],[0,0.201,-0.097],[1,0.216,-0.198],[1,0.211,-0.177],[1,0.217,-0.218],[2,-0.564,-0.206],[2,-0.561,-0.203],[2,-0.563,-0.207],[5,0.132,0.109],[5,0.145,0.09],[5,0.149,0.103],[3,0.223,-0.087],[3,0.215,-0.082],[3,0.226,-0.094],[6,-0.135,0.26],[6,-0.099,0.216],[6,-0.097,0.24],[4,-0.102,0.221],[4,-0.09,0.184],[4,-0.129,0.287]],"fp20":[[0,-0.066,-0.1],[0,-0.088,-0.28],[0,-0.076,-0.238],[1,-0.128,-0.252],[1,-0.139,-0.204],[1,-0.124,-0.204],[2,-0.128,0.206],[2,-0.104,0.028],[2,-0.084,0.2],[5,-0.104,0.106],[5,-0.132,0.065],[5,-0.165,0.124],[3,-0.148,-0.119],[3,-0.156,-0.056],[3,-0.136,-0.162],[6,0.012,0.366],[6,0.014,0.37],[6,0.01,0.34],[4,0.58,-0.061],[4,0.578,-0.062],[4,0.583,-0.068]],"fp40":[[0,-0.263,-0.097],[0,-0.305,-0.142],[0,-0.315,-0.157],[1,-0.136,-0.101],[1,-0.208,-0.096],[1,-0.203,-0.075],[2,-0.149,0.185],[2,-0.101,0.18],[2,-0.154,0.199],[5,0.303,-0.085],[5,0.318,-0.06],[5,0.31,-0.077],[3,0.226,-0.378],[3,0.225,-0.382],[3,0.221,-0.381],[6,0.27,0.377],[6,0.257,0.378],[6,0.257,0.375],[4,-0.186,0.11],[4,-0.179,0.113],[4,-0.19,0.114]]}

GROUPS = D["g"]
DEPTHS = ["0–20 cm","20–40 cm"]

# ---------- style ----------
plt.rcParams.update({
    "font.family": "Liberation Serif",
    "font.serif": ["Liberation Serif","Times New Roman","DejaVu Serif"],
    "mathtext.fontset": "stix",
    "font.size": 10.5,
    "axes.labelsize": 10.8,
    "xtick.labelsize": 8.3,
    "ytick.labelsize": 8.5,
    "axes.linewidth": 0.95,
    "xtick.direction": "out",
    "ytick.direction": "out",
    "xtick.major.width": 0.9,
    "ytick.major.width": 0.9,
    "xtick.major.size": 3.5,
    "ytick.major.size": 3.5,
})

BACT_COLORS = ['#4C78A8','#59A14F','#B279A2','#F28E2B','#EDC948','#76B7B2','#E15759','#9A9A9A']
FUNGI_COLORS = ['#4C78A8','#F28E2B','#59A14F','#B279A2','#EDC948','#76B7B2','#9A9A9A']
TREAT_COLORS = ['#1f77b4','#ff7f0e','#2ca02c','#d62728','#9467bd','#8c564b','#e377c2']

fig = plt.figure(figsize=(12.6, 15.6), facecolor="white")
outer = fig.add_gridspec(3,2,left=.035,right=.985,top=.985,bottom=.025,wspace=.055,hspace=.095)

def panel_frame(cell, letter, title):
    frame=fig.add_subplot(cell)
    frame.set_facecolor("white")
    frame.set_xticks([]); frame.set_yticks([])
    for sp in frame.spines.values():
        sp.set_visible(True); sp.set_linewidth(1.0); sp.set_color("black")
    frame.add_patch(Rectangle((0,.92),1,.08,transform=frame.transAxes,facecolor="#ededed",edgecolor="none",zorder=0))
    frame.plot([0,1],[.92,.92],transform=frame.transAxes,color="#777",lw=.8,clip_on=False)
    frame.text(.028,.96,f"({letter})",transform=frame.transAxes,ha="left",va="center",fontsize=13.0,fontweight="bold")
    frame.text(.50,.96,title,transform=frame.transAxes,ha="center",va="center",fontsize=13.2,fontweight="bold")
    return frame

def add_rel(frame, rect):
    b=frame.get_position()
    l,bo,w,h=rect
    return fig.add_axes([b.x0+b.width*l,b.y0+b.height*bo,b.width*w,b.height*h])

def close_box(ax):
    for s in ["left","right","top","bottom"]:
        ax.spines[s].set_visible(True); ax.spines[s].set_linewidth(.95); ax.spines[s].set_color("black")

def comp_panel(cell,letter,title,taxa,A,B,colors):
    frame=panel_frame(cell,letter,title)
    ax=add_rel(frame,[.11,.305,.865,.585])
    xl=np.arange(7); xr=np.arange(7)+8.25
    bl=np.zeros(7); br=np.zeros(7)
    hs=[]
    for i,tax in enumerate(taxa):
        vl=np.asarray(A[i],float); vr=np.asarray(B[i],float)
        h=ax.bar(xl,vl,bottom=bl,width=.74,color=colors[i],edgecolor="none")
        ax.bar(xr,vr,bottom=br,width=.74,color=colors[i],edgecolor="none")
        bl+=vl; br+=vr; hs.append(h[0])
    ax.set_ylim(0,100); ax.set_yticks(np.arange(0,101,20))
    ax.set_ylabel("Relative abundance (%)",fontweight="bold")
    ax.set_xticks(np.r_[xl,xr]); ax.set_xticklabels(GROUPS+GROUPS,rotation=45,ha="right",fontsize=7.9)
    ax.tick_params(axis="x",pad=.6)
    ax.axvline(7.62,color=".73",ls=(0,(3,3)),lw=.8)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)

    # depth bars sit close to rotated treatment labels
    ybar=-.155
    for xs,label in [(xl,"0–20 cm"),(xr,"20–40 cm")]:
        ax.plot([xs[0]-.30,xs[-1]+.30],[ybar,ybar],transform=ax.get_xaxis_transform(),color="black",lw=.8,clip_on=False)
        ax.text(xs.mean(),-.215,label,transform=ax.get_xaxis_transform(),ha="center",va="top",fontsize=9.8,fontweight="bold")

    # legend larger and pulled upward; no oversized empty band
    frame.legend(hs,taxa,ncol=4,loc="lower center",bbox_to_anchor=(.5,.018),frameon=False,
                 fontsize=8.6,handlelength=1.0,handletextpad=.35,columnspacing=.72,borderaxespad=0)

def shannon_panel(cell,letter,title,S,ymin,ymax,ystep):
    frame=panel_frame(cell,letter,title)
    ax=add_rel(frame,[.135,.165,.835,.715])
    x=np.arange(7); off=.16
    for dep,marker,col in [(0,"o",TREAT_COLORS[0]),(1,"^",TREAT_COLORS[1])]:
        means=np.array([S[dep][i][0] for i in range(7)])
        ses=np.array([S[dep][i][1] for i in range(7)])
        pos=x+(-off if dep==0 else off)
        ax.errorbar(pos,means,yerr=ses,fmt=marker,ls="none",ms=4.5,capsize=3.2,capthick=1.0,elinewidth=1.0,mew=.65,color=col,zorder=3)
    ax.set_xlim(-.25,6.25); ax.set_ylim(ymin,ymax); ax.set_yticks(np.arange(ymin,ymax+1e-9,ystep))
    ax.set_ylabel("Shannon index",fontweight="bold")
    ax.set_xticks(x); ax.set_xticklabels(GROUPS,rotation=45,ha="right",fontsize=8.2)
    ax.tick_params(axis="x",pad=1.0); close_box(ax)
    dh=[Line2D([0],[0],marker="o",linestyle="None",markersize=6,markerfacecolor=TREAT_COLORS[0],markeredgecolor=TREAT_COLORS[0],label="0–20 cm"),
        Line2D([0],[0],marker="^",linestyle="None",markersize=6.8,markerfacecolor=TREAT_COLORS[1],markeredgecolor=TREAT_COLORS[1],label="20–40 cm")]
    ax.legend(handles=dh,loc="upper left",frameon=False,fontsize=9.0,handletextpad=.5,borderaxespad=.6)

def pcoa_panel(cell,letter,title,left,right):
    frame=panel_frame(cell,letter,title)
    # enlarged square axes occupying more of panel
    rects=[[.10,.185,.385,.525],[.595,.185,.385,.525]]
    axes=[add_rel(frame,r) for r in rects]
    for ax,dep,dat,r in zip(axes,["0–20 cm","20–40 cm"],[left,right],rects):
        pts=np.asarray(dat["pts"],float)
        for i,g in enumerate(GROUPS):
            ss=pts[pts[:,0]==i]
            ax.scatter(ss[:,1],ss[:,2],s=27,color=TREAT_COLORS[i],edgecolor="black",linewidth=.30,zorder=3)
        xmin,xmax,xstep,ymin,ymax,ystep=dat["range"]
        ax.set_xlim(xmin,xmax); ax.set_xticks(np.arange(xmin,xmax+1e-9,xstep))
        ax.set_ylim(ymin,ymax); ax.set_yticks(np.arange(ymin,ymax+1e-9,ystep))
        ax.axhline(0,color=".74",ls=(0,(3,3)),lw=.75,zorder=0)
        ax.axvline(0,color=".74",ls=(0,(3,3)),lw=.75,zorder=0)
        ax.set_xlabel(f"PCoA1 ({dat['p1']:.1f}%)",fontsize=9.3,fontweight="bold",labelpad=1.2)
        ax.set_ylabel("")
        ax.tick_params(axis="both",labelsize=8.0,pad=1.2)
        close_box(ax)

        xc=r[0]+r[2]/2
        frame.add_patch(Rectangle((xc-.145,.812),.29,.047,transform=frame.transAxes,facecolor="#ededed",edgecolor="none"))
        frame.text(xc,.835,dep,transform=frame.transAxes,ha="center",va="center",fontsize=9.5,fontweight="bold")
        frame.text(xc,.770,f"Treatment: R$^2$={dat['r2']:.3f}, P=0.001",transform=frame.transAxes,ha="center",va="center",fontsize=8.2)

    # manual y labels: protected gutters inside outer panel
    frame.text(.048,.448,f"PCoA2 ({left['p2']:.1f}%)",rotation=90,transform=frame.transAxes,ha="center",va="center",fontsize=9.3,fontweight="bold")
    frame.text(.542,.448,f"PCoA2 ({right['p2']:.1f}%)",rotation=90,transform=frame.transAxes,ha="center",va="center",fontsize=9.3,fontweight="bold")

    hh=[Line2D([0],[0],marker="o",linestyle="None",markersize=5.8,markerfacecolor=TREAT_COLORS[i],markeredgecolor="black",markeredgewidth=.25,label=g) for i,g in enumerate(GROUPS)]
    frame.legend(handles=hh,ncol=7,loc="lower center",bbox_to_anchor=(.5,.030),frameon=False,
                 fontsize=8.0,handletextpad=.32,columnspacing=.62,borderaxespad=0,
                 title="Treatment",title_fontsize=8.6)

comp_panel(outer[0,0],"A","Bacterial dominant phyla",D["bt"],D["b20"],D["b40"],BACT_COLORS)
comp_panel(outer[0,1],"B","Fungal dominant phyla",D["ft"],D["f20"],D["f40"],FUNGI_COLORS)
shannon_panel(outer[1,0],"C","Bacterial Shannon",D["bs"],6.7,7.5,.1)
shannon_panel(outer[1,1],"D","Fungal Shannon",D["fs"],3.0,6.0,.5)

bp20={"pts":D["bp20"],"range":(-.4,.4,.2,-.3,.3,.1),"p1":15.0,"p2":9.8,"r2":.527}
bp40={"pts":D["bp40"],"range":(-.6,.4,.2,-.3,.3,.1),"p1":25.7,"p2":11.3,"r2":.648}
fp20={"pts":D["fp20"],"range":(-.2,.6,.2,-.4,.4,.2),"p1":20.7,"p2":14.2,"r2":.734}
fp40={"pts":D["fp40"],"range":(-.4,.4,.2,-.4,.4,.2),"p1":17.2,"p2":16.2,"r2":.798}
pcoa_panel(outer[2,0],"E","Bacterial PCoA",bp20,bp40)
pcoa_panel(outer[2,1],"F","Fungal PCoA",fp20,fp40)

out="fig3_cloud/output"
os.makedirs(out,exist_ok=True)
fig.savefig(f"{out}/Fig3_v12_3x2.png",dpi=450,facecolor="white",bbox_inches="tight",pad_inches=.02)
fig.savefig(f"{out}/Fig3_v12_3x2.pdf",facecolor="white",bbox_inches="tight",pad_inches=.02)
fig.savefig(f"{out}/Fig3_v12_3x2_600dpi.tif",dpi=600,facecolor="white",bbox_inches="tight",pad_inches=.02,
            pil_kwargs={"compression":"tiff_lzw"})
plt.close(fig)
print("done")
