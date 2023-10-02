#### Yu-Cheng 2022.09.04 modify
#### Jimmy 2022.05.25
#### Extract band data from vasp.
#### LNONCOLLINEAR = .TRUE. , LORBIT = 11
#### This is designed for noncollinear calculation.
#### Needed POSCAR for reciprocal conversion 
#### (that means only if band calculation is in reciprocal we can use this code)
#### 2020.05.07 Dealt with negative reciprocal vector by drecting reading blocks instaed of .split()
#### 2020.05.07 Output for xmgrace.
#### 2022.05.25 Rewritten.


import numpy as np
import math
import matplotlib.pyplot as plt

#### getting the basis vector for converting reciprocal
pos = open("/Users/yucheng/Desktop/Study/碩三上/Lab/data/CMA-c_GXUZX-GY/POSCAR","r")
print("Filename is :",pos.name)
print("Just eliminating:",pos.readline())
lattconst = float(pos.readline().split()[0])
a1,a2,a3 = [], [], []   #### a1-3 is basis vector from POSCAR
b1,b2,b3 = [], [], []
a1 = pos.readline().split()
a1 = [lattconst*float(x) for x in a1]
a2 = pos.readline().split()
a2 = [lattconst*float(x) for x in a2]
a3 = pos.readline().split()
a3 = [lattconst*float(x) for x in a3]
print(a1,a2,a3)
b1 = 2*math.pi*np.cross(a2,a3)/np.inner(a1,np.cross(a2,a3))   
b2 = 2*math.pi*np.cross(a3,a1)/np.inner(a2,np.cross(a3,a1))
b3 = 2*math.pi*np.cross(a1,a2)/np.inner(a3,np.cross(a1,a2))
print(b1,b2,b3)
pos.close()

#### Wished atoms and dos orbitals
#### input format
# 1. What atoms do you need? (start from 1 in POSCAR order)
# 2. Projection? py=1 pz=2 px=3 dxy=4 dyz=5 dz2=6 dxz=7 dx2=8 
#                m1(dxz, dyz)=9 m2(dxy, dx2)=10, eg=11, t2g=12
#  (s,p,d will be out by default.)
# 3. output name
b_input = open("/Users/yucheng/Desktop/Study/碩三上/Lab/data/CMA-c_GXUZX-GY/Band_input.txt")
iatom = b_input.readline().split()
iatom = [int(x) for x in iatom ]
iproj = b_input.readline().split()
iproj = [int(x) for x in iproj ]
outname = b_input.readline()
b_input.close()

##### ask for wanted bands
#efermi = float(input("Fermi energy?\n"))
efermi = 5.71442306 # from vasprun.xml of 0_scf calculation; LAG: 7.15101731, CMA-c: 5.71442306
#print("There are ", nband,"bands.\n")
bmin = int(input("First band?(band index starting from 1)\n"))
bmax = int (input("Last band?\n"))
wband = bmax-bmin+1 #for wanted bands
wbandindex = []
for i in range(wband):
    wbandindex.append(i+bmin)

print(wbandindex)

#### getting number of ions, number of kpoints, number of bands 
band = open("/Users/yucheng/Desktop/Study/碩三上/Lab/data/CMA-c_GXUZX-GY/PROCAR","r")
print("Filename is :",band.name)
print("Just eliminating",band.readline())
tmp = band.readline()
nkpt = int(tmp.split()[3])
nband = int(tmp.split()[7])
natom = int(tmp.split()[11])

### Confirmation of all used data
print("efermi:",efermi,"nkpoints:",nkpt,"nbands:",nband,"wanted bands:",wband,"natoms;",natom)
print("band index:",wbandindex)
print("Filename will be in",outname)

#### Matrix definition  one array for one band, and nkpt slots plus 4 [Tot, Sx, Sy, Sz]                                           
kpt,energy = np.zeros((wband,nkpt)),np.zeros((wband,nkpt))
s = np.zeros((wband,nkpt,4)) 
py,pz,px = np.zeros((wband,nkpt,4)),np.zeros((wband,nkpt,4)),np.zeros((wband,nkpt,4))
p = np.zeros((wband,nkpt,4))
dxy,dx2 = np.zeros((wband,nkpt,4)),np.zeros((wband,nkpt,4))
dyz,dxz = np.zeros((wband,nkpt,4)),np.zeros((wband,nkpt,4))
dz2 = np.zeros((wband,nkpt,4))  # m=0
dm1 = np.zeros((wband,nkpt,4)) # m=1 , dyz+dxz          
dm2 = np.zeros((wband,nkpt,4)) # m=2 , dxy,dx2-y2
eg, t2g = np.zeros((wband,nkpt,4)),np.zeros((wband,nkpt,4))
d = np.zeros((wband,nkpt,4))
all_atom = np.zeros((wband,nkpt,4))
highsympot = []

print("Start reading")
#### kpoint definition : a = a + distance between one and the one before
#### kpoint matrix fill in 
kdis = 0.0         # kpoint distance sum
kkdis = 0.0        # distance between k_now and prev_k
nowk = np.array([0.0,0.0,0.0])
prevk = np.array([0.0,0.0,0.0]) 
band.readline() #line 3
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
    elif tmp.split()[0]=="k-point":
        k_now = int(tmp[8:14])-1       # kpoint number k_now
        k1 = float(tmp[18:29])         # Direct coordinate k1, k2, k3 
        k2 = float(tmp[29:40])
        k3 = float(tmp[40:52])
        nowk = k1*b1+k2*b2+k3*b3          # Change to cartesian by timing b1,b2,b3
        kkdis = math.sqrt((nowk[0]-prevk[0])**2+(nowk[1]-prevk[1])**2+(nowk[2]-prevk[2])**2)
        kdis = kdis + kkdis                  # 和當前k(nowk)和前一個k(prevk)的距離 再一一加上去
        for i in range(wband):
            kpt[i][k_now] = kdis 
        if kkdis == 0:
            highsympot.append(kdis)     # since in KPOINTS we repear high symmetry points, it also repeat in PROCAR
        prevk = nowk
    elif tmp.split()[0]=="band":             #直接繼續往下讀
        b_now = int(tmp.split()[1])          # what band we are on
        if b_now in wbandindex:
            energy[b_now-1][k_now]= float(tmp.split()[4])-efermi  # since energy has shape=(wband,nkpt), you can only start
                                                                  # from band index=1 or index out of range
            band.readline()
            band.readline()
            for k in range(4):
                for j in range(natom):     #可以多讀nbands行 因為我們確定在band number後前面是全部貢獻 後面才是 mx my mz
                    read_atom = band.readline()
                    # print("check-point-1", "b_now= ", b_now, "k= ", k, "j= ", j, "read_atom.split()= ", read_atom.split())
                    if int(read_atom.split()[0]) in iatom:
                        s[b_now-1][k_now][k] = s[b_now-1][k_now][k]+float(read_atom.split()[1])             
                        py[b_now-1][k_now][k] = py[b_now-1][k_now][k]+float(read_atom.split()[2])           
                        pz[b_now-1][k_now][k] = pz[b_now-1][k_now][k]+float(read_atom.split()[3])           
                        px[b_now-1][k_now][k] = px[b_now-1][k_now][k]+float(read_atom.split()[4])           
                        p[b_now-1][k_now][k]=py[b_now-1][k_now][k]+pz[b_now-1][k_now][k]+px[b_now-1][k_now][k]  
                        dxy[b_now-1][k_now][k] = dxy[b_now-1][k_now][k]+float(read_atom.split()[5])        
                        dyz[b_now-1][k_now][k] = dyz[b_now-1][k_now][k]+float(read_atom.split()[6])         
                        dz2[b_now-1][k_now][k] = dz2[b_now-1][k_now][k]+float(read_atom.split()[7])       
                        dxz[b_now-1][k_now][k] = dxz[b_now-1][k_now][k]+float(read_atom.split()[8])
                        dx2[b_now-1][k_now][k] = dx2[b_now-1][k_now][k]+float(read_atom.split()[9])
                        dm1[b_now-1][k_now][k] = dyz[b_now-1][k_now][k]+dxz[b_now-1][k_now][k]
                        dm2[b_now-1][k_now][k] = dxy[b_now-1][k_now][k]+dx2[b_now-1][k_now][k]
                        eg[b_now-1][k_now][k] = dz2[b_now-1][k_now][k]+dx2[b_now-1][k_now][k]
                        t2g[b_now-1][k_now][k] = dxy[b_now-1][k_now][k]+dyz[b_now-1][k_now][k]+dxz[b_now-1][k_now][k]
                        d[b_now-1][k_now][k]=dxy[b_now-1][k_now][k]+dyz[b_now-1][k_now][k]+dz2[b_now-1][k_now][k]+dxz[b_now-1][k_now][k]+dx2[b_now-1][k_now][k]
                read_atom = band.readline()
                # print("check-point-2", "read_atom= ", read_atom)
                all_atom[b_now-1][k_now][k] = float(read_atom.split()[10])
                # band.readline()
band.close()

std_out = {0:"s",1:"p",2:"d"}
std_out_dat = {0:s,1:p,2:d}
#### Dictionary to store name of orbitals
orbitals = {1:"py",2:"pz",3:"px",4:"dxy",5:"dyz",6:"dz2",7:"dxz",8:"dx2",9:"dm1",10:"dm2",11:"eg",12:"t2g"}
#### Dictionary that contains arrays for projected orbitals  
orb_dat = {1:py,2:pz,3:px,4:dxy,5:dyz,6:dz2,7:dxz,8:dx2,9:dm1,10:dm2,11:eg,12:t2g}
#### Dictionary containing spin rpojection direction
sp_direc = {0:"Tot", 1:"Sx", 2:"Sy", 3:"Sz"}

sp_proj = int(input("Spin projection? 0-3 : Tot Sx Sy Sz \n"))  # Tot, Sx, Sy, Sz         
#### Write the needed files out
out = {}
out['band'] = open(outname+"-band.dat","w")  # why don't just open 
for i in range(3):
    out[std_out[i]] =open(outname+"-"+std_out[i]+"-"+sp_direc[sp_proj]+".dat","w")
# First output high symmetry position
for i in range(len(highsympot)):
    out['band'].write(str(highsympot[i])+" ")
    for j in range(3):
        out[std_out[j]].write(str(highsympot[i])+" ")
out['band'].write("\n")
for i in range(3):
    out[std_out[i]].write("\n")
# Second output the band information
for i in range(wband):
    for j in range(nkpt):
        out['band'].write(str(kpt[i][j])+" "+str(energy[i][j])+"\n")
        for k in range(3):
            out[std_out[k]].write(str(kpt[i][j])+" "+str(energy[i][j])+" "+str(std_out_dat[k][i][j][sp_proj])+"\n")
#Seperate each bands    
    out['band'].write("\n")
    for l in range(3):
        out[std_out[l]].write("\n") 
out['band'].close()
for m in range(3):
        out[std_out[m]].close()

if 0 not in iproj:
    for i in iproj:
        out[orbitals[i]] =open(outname+"-"+orbitals[i]+"-"+sp_direc[sp_proj]+".dat","w")
    for i in iproj:
        for j in range(len(highsympot)):
            out[orbitals[i]].write(str(highsympot[j])+" ")
        out[orbitals[i]].write("\n")    
    for i in range(wband):
        for j in range(nkpt):
            for k in iproj:
                out[orbitals[k]].write(str(kpt[i][j])+" "+str(energy[i][j])+" "+str(orb_dat[k][i][j][sp_proj])+"\n")
        for l in iproj:
            out[orbitals[l]].write("\n")
    for m in iproj :
        out[orbitals[m]].close()

print("High symmetry at: ",highsympot)   
print("All finish!")
