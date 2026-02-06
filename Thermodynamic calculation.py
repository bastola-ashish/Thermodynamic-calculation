# %%

# =============================================================================
# Import required packages
# =============================================================================

from reaktoro import *
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# =============================================================================
# Import results saved in a file
# =============================================================================

icp = pd.read_excel("result.xlsx",sheet_name="icp")    # ICP result
strength = pd.read_excel("result.xlsx",sheet_name="strength")   # UCS result

# =============================================================================
# Database selection
# Ensure that the database file is in the working directory
# =============================================================================

db = ThermoFunDatabase("cemdata18")

# =============================================================================
# System definition and activity model selection
# =============================================================================

aqueous = AqueousPhase(speciate("Al Si Ca K Mg Na S O H"))

system = ChemicalSystem(db, aqueous.set(ActivityModelPitzer()))

# =============================================================================
# Define specification to setup pH constraint
# =============================================================================

# =============================================================================
# Speciation Calculation
# =============================================================================
chemicalprop_list=[]
for i in range(len(icp)):
    specs = EquilibriumSpecs(system)
    specs.temperature()
    specs.pressure()
    specs.pH()

    conditions = EquilibriumConditions(specs)
    conditions.temperature(25.0, "celsius")
    conditions.pressure(1.0, "atm")
    conditions.pH(icp['pH'][i])

    state = ChemicalState(system)
    state.set("H2O@",  1, "kg")
    state.set("Al+3",  icp['Al'][i], "mol")
    state.set("Ca+2",  icp['Ca'][i], "mol")
    state.set("K+",  icp['K'][i], "mol")
    state.set("Na+",  icp['Na'][i], "mol")
    state.set("SO4-2",  icp['S'][i], "mol")
    state.set("SiO2@",  icp['Si'][i], "mol")

    # Define solver
    solver = EquilibriumSolver(specs)

    # Solve chemical state for 0 min
    solver.solve(state, conditions)
    chemicalprop_list.append(ChemicalProps(state))
    

# %%
# =============================================================================
# Saturation index of the phases of interest
# =============================================================================
icp["csh"]=[max(aprops.saturationIndex("CSH3T-T2C"),aprops.saturationIndex("CSH3T-T5C"),
        aprops.saturationIndex("CSH3T-TobH"),aprops.saturationIndex("CSHQ-JenD"),
        aprops.saturationIndex("CSHQ-JenH"),aprops.saturationIndex("CSHQ-TobD"),
        aprops.saturationIndex("CSHQ-TobH"))[0] for aprops in 
        [AqueousProps(chemicalprop) for chemicalprop in chemicalprop_list]]

icp["ettrin"] = [max(aprops.saturationIndex("ettringite"),aprops.saturationIndex("ettringite13"),
              aprops.saturationIndex("Ettringite13_des"),aprops.saturationIndex("ettringite30"),
              aprops.saturationIndex("ettringite9"),aprops.saturationIndex("Ettringite9_des"))[0] for aprops in 
              [AqueousProps(chemicalprop) for chemicalprop in chemicalprop_list]]

icp["gyp"]=[aprops.saturationIndex("Gp")[0] for aprops in 
        [AqueousProps(chemicalprop) for chemicalprop in chemicalprop_list]]

icp["gib"]=[aprops.saturationIndex("Gbs")[0] for aprops in 
        [AqueousProps(chemicalprop) for chemicalprop in chemicalprop_list]]

icp["port"]=[aprops.saturationIndex("Portlandite")[0] for aprops in 
        [AqueousProps(chemicalprop) for chemicalprop in chemicalprop_list]]

icp["strat"]=[max(aprops.saturationIndex("straetlingite"),aprops.saturationIndex("straetlingite5.5"),
          aprops.saturationIndex("straetlingite7"))[0] for aprops in 
          [AqueousProps(chemicalprop) for chemicalprop in chemicalprop_list]]

icp["mono"]=[max(aprops.saturationIndex("monosulphate10.5"),aprops.saturationIndex("monosulphate12"),
          aprops.saturationIndex("monosulphate1205"),aprops.saturationIndex("monosulphate14"),
          aprops.saturationIndex("monosulphate16"),aprops.saturationIndex("monosulphate9"))[0] for aprops in 
          [AqueousProps(chemicalprop) for chemicalprop in chemicalprop_list]]

icp["csh_eff"]=[max(aprops.saturationIndex("CSH3T-T2C")/5.5,aprops.saturationIndex("CSH3T-T5C")/5,
        aprops.saturationIndex("CSH3T-TobH")/4.5,aprops.saturationIndex("CSHQ-JenD")/5.167,
        aprops.saturationIndex("CSHQ-JenH")/4.999,aprops.saturationIndex("CSHQ-TobD")/3.166825,
        aprops.saturationIndex("CSHQ-TobH")/3.0001)[0] for aprops in 
        [AqueousProps(chemicalprop) for chemicalprop in chemicalprop_list]]

icp["gyp_eff"]=[aprops.saturationIndex("Gp")[0]/2 for aprops in 
        [AqueousProps(chemicalprop) for chemicalprop in chemicalprop_list]]

icp["gib_eff"]=[aprops.saturationIndex("Gbs")[0]/2 for aprops in 
        [AqueousProps(chemicalprop) for chemicalprop in chemicalprop_list]]

icp["port_eff"]=[aprops.saturationIndex("Portlandite")[0]/3 for aprops in 
        [AqueousProps(chemicalprop) for chemicalprop in chemicalprop_list]]

icp["strat_eff"]=[max(aprops.saturationIndex("straetlingite"),aprops.saturationIndex("straetlingite5.5"),
          aprops.saturationIndex("straetlingite7"))[0]/7 for aprops in 
          [AqueousProps(chemicalprop) for chemicalprop in chemicalprop_list]]

icp["ettrin_eff"] = [max(aprops.saturationIndex("ettringite"),aprops.saturationIndex("ettringite13"),
              aprops.saturationIndex("Ettringite13_des"),aprops.saturationIndex("ettringite30"),
              aprops.saturationIndex("ettringite9"),aprops.saturationIndex("Ettringite9_des"))[0]/15 for aprops in 
              [AqueousProps(chemicalprop) for chemicalprop in chemicalprop_list]]

icp["mono_eff"]=[max(aprops.saturationIndex("monosulphate10.5"),aprops.saturationIndex("monosulphate12"),
          aprops.saturationIndex("monosulphate1205"),aprops.saturationIndex("monosulphate14"),
          aprops.saturationIndex("monosulphate16"),aprops.saturationIndex("monosulphate9"))[0]/11 for aprops in 
          [AqueousProps(chemicalprop) for chemicalprop in chemicalprop_list]]

# %%
# =============================================================================
# Plot the result
# =============================================================================

fig,axes=plt.subplots(2,2,figsize=[6.5,6.5])
axes_twin = [ax.twinx() for ax in axes.flatten()]
axes_twin_second = [ax.twinx() for ax in axes.flatten()]

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial"]
plt.rcParams["font.size"]=12

for ax in axes_twin:
    ax.set_ylim([0,300])
    ax.spines.top.set_visible(0)
    ax.spines.left.set_visible(0)
    ax.spines.bottom.set_visible(0)
    ax.spines.right.set_color("red")
    ax.tick_params(axis='y', colors='red')
    ax.set_yticks(np.linspace(0,300,5))
    if ax in [axes_twin[0],axes_twin[2]]:
        ax.set_yticks(np.linspace(0,300,5),[])

for ax in axes_twin_second:
    ax.set_ylim([4,12])
    ax.set_yticks(np.linspace(4,12,5))
    ax.spines.top.set_visible(0)
    ax.spines.left.set_visible(0)
    ax.spines.bottom.set_visible(0)
    ax.spines["right"].set_position(("axes", 1.22))
    ax.spines["right"].set_color("blue")
    ax.tick_params(axis='y', colors='blue')
    if ax in [axes_twin_second[0],axes_twin_second[2]]:
        ax.spines["right"].set_position(("axes",1.1))
        ax.set_yticks(np.linspace(4,12,5),[])
        # ax.spines.right.set_visible(0)
        # ax.set_yticks([],[])


for ax in axes.flatten():
    ax.grid(lw=0.75,alpha=0.5)
    ax.set_xlim([0,6])
    ax.hlines(0,6,0,color="r",ls="--")
    ax.spines.right.set_visible(0)
    ax.set_ylim([-6,2])
    
    
x=np.arange(6)
x1 = np.array([0,3,4,5])

ax = axes[0,1]
ax_1 = axes_twin[1]
data = icp[icp["Stabilizer"]=="csa"]
ax.plot(x,data["csh_eff"],"k--o",mfc="w",label="C-S-H")
ax.plot(x,data["ettrin_eff"],"k--d",mfc="w",label="Ettringite")
ax.plot(x,data["strat_eff"],"k--x",mfc="w",label="Straetlingite")
ax.plot(x,data["gyp_eff"],"k--v",mfc="w",alpha=0.5,label="Gypsum")
ax.plot(x,data["gib_eff"],"k--^",mfc="w",alpha=0.5,label="Gibbsite")
ax.plot(x,data["port_eff"],"k--o",mfc="w",alpha=0.5,label="Portlandite")
ax.plot(x,data["mono_eff"],"k--",mfc="w",alpha=0.5,label="Monosulfate")
ax_1.plot(x1,strength["csa"]*100,"r--o",mfc="w")
h1,l1 = ax.get_legend_handles_labels()
h2,l2 = ax_1.get_legend_handles_labels()
ax.legend(handles=h1+h2,labels=l1+l2,ncols=2,fontsize=8)


plt.tight_layout()

# %%



x=np.arange(1,7,1)

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial"]
plt.rcParams["font.size"]=10
plt.rcParams["font.style"]="italic"
plt.rcParams["font.weight"]="bold"
plt.rcParams["lines.linewidth"]=0.75

# define axis behavior

fig , ax1 = plt.subplots(figsize=[4.7,3])
ax2 = ax1.twinx()
ax2.spines["right"].set_color("red")
ax2.tick_params(axis='y', colors='red') 
ax3 = ax1.twinx()
ax3.spines["right"].set_position(("axes", 1.22))
ax3.spines["right"].set_color("blue")
ax3.tick_params(axis='y', colors='blue') 

# plot result

line1=ax1.plot(x,csh_eff,'ko--',markersize=3,fillstyle='none',markeredgewidth=0.5,label="C-S-H",mfc='k')
line2=ax1.plot(x,ettrin_eff,'kd--',markersize=3,fillstyle='none',markeredgewidth=0.5,label="Ettringite",mfc='k')
line3=ax1.plot(x,strat_eff,'kx--',markersize=3,label="Stratlingite")

line4=ax1.plot(x,gyp_eff,'o--',color='0.5',markersize=3,label="Gypsum")
line5=ax1.plot(x,gib_eff,'v--',color='0.5',markersize=3,label="Gibbsite")
line6=ax1.plot(x,port_eff,'^--',color='0.5',markersize=3,label="Portlandite")
line7=ax1.plot(x,mono_eff,'x-',color='0.5',markersize=3,label="Monosulfate")

ax1.set_xlabel("Curing duration",weight="bold")
ax1.set_xlim([0.5,6.5])
ax1.set_xticks(np.arange(1,7,1),labels=np.array(['0 min','15 min','30 min','1 h','1 d','7 d']))
ax1.set_ylim([-6,2])
ax1.set_ylabel("Effective saturation index",weight="bold")
# ax1.legend(ncol=2,fontsize=8,loc="lower right",frameon=False,borderpad=0.1)

ax1.grid(lw=0.1,alpha=.7)

x2 = np.array([1,4,5,6])
line8=ax2.plot(x2,strength*100,'r--.',label="Strength",markeredgewidth=0.5,fillstyle='none')
ax2.set_ylim([0,300])
ax2.set_ylabel("Relative UCS (%)",weight="bold",color='r')
# ax2.legend(ncol=2,fontsize=8,loc="lower right",frameon=False,borderpad=0.1)

line9=ax3.plot(x,df['pH'],'b--.',fillstyle='none',markeredgewidth=0.5,label="pH")
ax3.set_ylabel("pH",weight="bold",color='b')
ax3.set_ylim([4,12])
# ax3.legend(ncol=2,fontsize=8,loc="right",frameon=False,borderpad=0.1)


line10=ax1.plot([4],[ettrin_eff[3]],'rd',markersize=3,label="Observed")

line = line1 + line2 + line3 + line4 + line5+ line6 + line7 + line8 + line9 + line10
labels = [l.get_label() for l in line]
ax1.legend(line, labels, ncol=2,fontsize=8,loc="lower right",frameon=False,borderpad=0.1)
# ax1.grid(lw=0.1,alpha=0.8)
plt.tight_layout()

# plt.savefig("stab_b.svg")

# %%
