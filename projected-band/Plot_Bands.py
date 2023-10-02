# 2022.09.04 Modify by Yu-Cheng
# 2022.07.11 Read from our own output and plot projected bands
# Input format
# L1 : High symmetry point
# k E-Ef proj *one band / bands info are separated with a space

# /Users/yucheng/Desktop/Study/碩三暑假/Lab/code/jimmy_pband/CuMnAs-soc_Mn-dz2-Sz.dat
import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from numpy.lib.function_base import append

## Interaction:
#efermi = input("Efermi?\n")
#bmin = input("Band min index? \n")
#bmax = input("Band max index? \n")
#openfile = ("File name? \n")

## Start reading PROCAR
efermi = 5.71442306  # LAG: 7.15101731, CMA-c: 5.71442306, 
bmin = 1
bmax = 120
wband = bmax-bmin+1 #for wanted bands
wbandindex = []
for i in range(wband):
    wbandindex.append(i+bmin)   # This is band index!!

high_symm_label=[r'$\Gamma$',"X","U","Z","X"]  # CMA-c
# high_symm_label=[r'$\Gamma$',r'$\Sigma$',"N",r'$\Sigma_1$',"Z",r'$\Gamma$',"X"]  # LAG
# plt.title(r'$\Gamma$')
# plt.show()

#### Matrix definition
kpt = []
energy = []
proj = []

openfile = str(input("What file? \n"))
band = open(openfile,"r")
print("Filename is :",band.name)
# The first line is the high symetry point
highsympot = band.readline().split()
while True:
    tmp = band.readline()
    if not tmp:  #check if out of file or reach EOF(end of file)
        break
    elif tmp == "\n":    #eliminating blanks due to error of list out of range for split(0)
        #print("blank")
        continue
    elif tmp == " \n":
        #print("blank")
        continue
    else:  # k E-Ef proj
        kpt.append(float(tmp.split()[0]))
        energy.append(float(tmp.split()[1]))
        proj.append(float(tmp.split()[2]))
band.close()

nkpt=int(len(kpt)/wband)  # since kpt is for all bands

tnrfont = {'fontname':'Times New Roman'}

print("Start plot!")
## Concept to draw band structure is that you create discrete data set and  you do a line segment collection and draw all the lines out.
## The projected method is then an average between the data and the k point next to it. 
x = np.zeros(nkpt)
y = np.zeros(nkpt)
proj_dat = np.zeros(nkpt-1)
proj_plot = np.empty([1])  # why we want to designate initial non-zero value? it will be covered up anyway
bigsegments =np.empty([2,2,2])
print(min(proj),max(proj),"\n")
#cbmin = input("Color bar min?\n")
#cbmax = input("Color bar max?\n")
#! cbmin = 0
#! cbmax = 0.01
# figure size
#plt.figure(figsize=(6.6,5.7))
fig, axs = plt.subplots(figsize=(9,6))
# read from each data
for i in range(wband):
    for j in range(nkpt):
        x[j]=kpt[j+i*nkpt]
        y[j]=energy[j+i*nkpt]
    # Projection average    
    for k in range(nkpt-1):    
        proj_dat[k] = (float(proj[k+i*nkpt])+float(proj[k+1+i*nkpt]))/2
        #print("sp:",surf_proj)
    proj_plot=np.append(proj_plot,proj,axis=0)  # why don't just use proj? --> extra dimension
    # Reformation into [ [[k1,E1]], [[k2,E2]], [[k3,E3]]  ]
    points = np.array([x, y]).T.reshape(-1, 1, 2)  # -1 means the len of that dimension is set according to the rest
    # Then every two terms are combined [ [[k1,E1],[k2,E2]], [[k2,E2],[k3,E3]] ]
    segments = np.concatenate([points[:-1], points[1:]], axis=1)  # [:-1] means from start to N-1
    # Append line collention of one band into big segment  
    bigsegments=np.append(bigsegments,segments,axis=0)
# Cancel the default first two
#bigsegments=np.delete(bigsegments,0,0)
#bigsegments=np.delete(bigsegments,0,0)

# Add projected data onto line collection 
proj_plot=np.delete(proj_plot,0)  # didn't use averaged projection value (proj_dat), instead remove first value
                                  # should be approximately euqal if k-points large enough

# Set the values used for colormapping
#! norm = plt.Normalize(cbmin,cbmax)  # this is not used as well
lc = LineCollection(bigsegments, colors="k", alpha =1.0, linewidths=2*proj_plot)  # why don't set linewidths=proj_plot here?

line = axs.add_collection(lc)

# Color bar settings 
#cb = plt.colorbar(line, ax=axs)
#cb.set_label(label='Weights', size=18,loc='center',rotation='vertical',**tnrfont)
#cb.ax.tick_params(labelsize=16)
#for l in cb.ax.yaxis.get_ticklabels():
#    l.set_family("Times New Roman")

# Plot high symmetry line out
axs.plot((x.min(),x.max()),(0,0),color='brown',linewidth=1,linestyle='--')  # (x_start,x_end), (y_start,y_end)
for i in range(len(highsympot)):
    axs.plot((float(highsympot[i]),float(highsympot[i])),(-12,5),color='brown',linewidth=1,linestyle='--') 

# High symmetry label
labelposxshift = 0.05
ylimlow = -1.5
ylimhigh = 1.5
labelpos = ylimlow-0.22
#from matplotlib import rc
#rc('font', **{'family':'sans-serif','sans-serif':['Times New Roman']})
#rc('text', usetex=True)
axs.text(float(min(kpt))-labelposxshift, labelpos, high_symm_label[0], size = 22, color = 'k',fontfamily='Times New Roman')
for i in range(len(high_symm_label)-2):
    axs.text(float(highsympot[i+1])-labelposxshift,labelpos, high_symm_label[i+1], size = 22, color = 'k',fontfamily='Times New Roman')
axs.text(float(max(kpt))-labelposxshift, labelpos, high_symm_label[-1], size = 22, color = 'k',fontfamily='Times New Roman')

# General plot and axis setteings
plt.tick_params(
axis='x',          # changes apply to the x-axis
which='major',      # both major and minor ticks are affected
bottom=False,      # ticks along the bottom edge are off
top=False,         # ticks along the top edge are off
labelbottom=False,
labelsize=18) # labels along the bottom edge are off
plt.tick_params(axis='y',which='major',labelsize=22)
for tick in axs.get_yticklabels():
    tick.set_fontname("Times New Roman")
#plt.yticks()
axs.set_xlabel(None)
axs.set_ylabel('Energy ( eV )',fontsize=22,**tnrfont)
axs.set_xlim(x.min(), x.max())                   ## Define upper and lower limit of x and y 
axs.set_ylim(ylimlow,ylimhigh)
#cb = plt.colorbar(fig)
#plt.title(title)
plt.savefig('Mn-dz2-Sz-orbital.png', dpi=1200)
# plt.show()    

print("All finish!")
