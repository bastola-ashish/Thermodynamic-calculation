# %%

# =============================================================================
# Import required packages
# =============================================================================

from reaktoro import *
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({
    "text.usetex": True,
    "text.latex.preamble": r"""
        \usepackage{helvet}       % Helvetica font
        \renewcommand{\familydefault}{\sfdefault} % Sans-serif default
        \usepackage{amsmath}     % AMS math
        \usepackage{sfmath}      % Sans-serif math
        \usepackage{bm}          % Bold math symbols
        \boldmath               
    """,
    "font.family": "sans-serif",
    "font.sans-serif": "Helvetica",
    "font.size": 10,
})

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

fig,axes=plt.subplots(2,2,figsize=[7,6.5])

axes_twin = [ax.twinx() for ax in axes.flatten()]
axes_twin_second = [ax.twinx() for ax in axes.flatten()]

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
    ax.set_xlim([0,5])
    ax.set_xticks(np.linspace(0,5,6),
                  ["0 min","15 min","30 min","1h","1d","7d"])
    ax.hlines(0,6,0,color="g",ls="--")
    ax.spines.right.set_visible(0)
    ax.set_ylim([-6,2])
    
    
x=np.arange(6)
x1 = np.array([0,3,4,5])

ax = axes[0,0]
ax_1 = axes_twin[0]
ax_2 = axes_twin_second[0]
data = icp[icp["Stabilizer"]=="type3"]
ax.plot(x,data["csh_eff"],"k--o",mfc="w",ms=4,label="C-S-H")
ax.plot(x,data["ettrin_eff"],"r--d",mfc="w",ms=4,label="Ettringite")
ax.plot(x,data["strat_eff"],"k-",mfc="w",ms=4,label="Straetlingite")
ax.plot(x,data["gyp_eff"],"k--v",mfc="w",ms=4,alpha=0.5,label="Gypsum")
ax.plot(x,data["gib_eff"],"k--^",mfc="w",ms=4,alpha=0.5,label="Gibbsite")
ax.plot(x,data["port_eff"],"k--o",mfc="w",ms=4,alpha=0.5,label="Portlandite")
ax.plot(x,data["mono_eff"],"k--",mfc="w",ms=4,alpha=0.5,label="Monosulfate")
ax_1.plot(x1,strength["type_3"]*100,"r--o",ms=4,mfc="w")
ax_2.plot(x,data["pH"],"b--o",ms=4,mfc="w")
ax.set_ylabel(r"\textbf{Effective Saturation Index}")

ax = axes[0,1]
ax_1 = axes_twin[1]
ax_2 = axes_twin_second[1]
data = icp[icp["Stabilizer"]=="csa"]
ax.plot(x,data["csh_eff"],"k--o",mfc="w",ms=4,label="C-S-H")
ax.plot(x,data["ettrin_eff"],"r--d",mfc="w",ms=4,label="Ettringite")
ax.plot(x,data["strat_eff"],"k-",mfc="w",ms=4,label="Straetlingite")
ax.plot(x,data["gyp_eff"],"k--v",mfc="w",ms=4,alpha=0.5,label="Gypsum")
ax.plot(x,data["gib_eff"],"k--^",mfc="w",ms=4,alpha=0.5,label="Gibbsite")
ax.plot(x,data["port_eff"],"k--o",mfc="w",ms=4,alpha=0.5,label="Portlandite")
ax.plot(x,data["mono_eff"],"k--",mfc="w",ms=4,alpha=0.5,label="Monosulfate")
ax_1.plot(x1,strength["type_3"]*100,"r--o",ms=4,mfc="w",label="UCS")
ax_2.plot(x,data["pH"],"b--o",ms=4,mfc="w",label="pH")
h1,l1 = ax.get_legend_handles_labels()
h2,l2 = ax_1.get_legend_handles_labels()
h3,l3 = ax_2.get_legend_handles_labels()
ax.legend(handles=h1+h2+h3,labels=l1+l2+l3,ncols=2,fontsize=8)
ax_1.set_ylabel(r"\textbf{UCS}",color="r")
ax_2.set_ylabel(r"\textbf{pH}",color="b")

ax = axes[1,0]
ax_1 = axes_twin[2]
ax_2 = axes_twin_second[2]
data = icp[icp["Stabilizer"]=="StabA"]
ax.plot(x,data["csh_eff"],"k--o",mfc="w",label="C-S-H")
ax.plot(x,data["ettrin_eff"],"k--d",mfc="w",label="Ettringite")
ax.plot(x,data["strat_eff"],"k--x",mfc="w",label="Straetlingite")
ax.plot(x,data["gyp_eff"],"k--v",mfc="w",alpha=0.5,label="Gypsum")
ax.plot(x,data["gib_eff"],"k--^",mfc="w",alpha=0.5,label="Gibbsite")
ax.plot(x,data["port_eff"],"k--o",mfc="w",alpha=0.5,label="Portlandite")
ax.plot(x,data["mono_eff"],"k--",mfc="w",alpha=0.5,label="Monosulfate")
ax_1.plot(x1,strength["stab_a"]*100,"r--o",mfc="w")
h1,l1 = ax.get_legend_handles_labels()
h2,l2 = ax_1.get_legend_handles_labels()
ax.set_ylabel(r"\textbf{Effective Saturation Index}")
ax.set_xlabel(r"\textbf{Curing duration}")

ax = axes[1,1]
ax_1 = axes_twin[3]
ax_2 = axes_twin_second[3]
data = icp[icp["Stabilizer"]=="StabB"]
ax.plot(x,data["csh_eff"],"k--o",mfc="w",label="C-S-H")
ax.plot(x,data["ettrin_eff"],"k--d",mfc="w",label="Ettringite")
ax.plot(x,data["strat_eff"],"k--x",mfc="w",label="Straetlingite")
ax.plot(x,data["gyp_eff"],"k--v",mfc="w",alpha=0.5,label="Gypsum")
ax.plot(x,data["gib_eff"],"k--^",mfc="w",alpha=0.5,label="Gibbsite")
ax.plot(x,data["port_eff"],"k--o",mfc="w",alpha=0.5,label="Portlandite")
ax.plot(x,data["mono_eff"],"k--",mfc="w",alpha=0.5,label="Monosulfate")
ax_1.plot(x1,strength["stab_b"]*100,"r--o",mfc="w")
h1,l1 = ax.get_legend_handles_labels()
h2,l2 = ax_1.get_legend_handles_labels()
ax_1.set_ylabel(r"\textbf{UCS}",color="r")
ax_2.set_ylabel(r"\textbf{pH}",color="b")
ax.set_xlabel(r"\textbf{Curing duration}")

plt.tight_layout()