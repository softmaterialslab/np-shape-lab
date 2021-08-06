import numpy as np
import csv
import scipy.linalg as la

# define spherical harmonics constants
Y00 = 0.5*(1/np.pi)**0.5
Y01 = lambda t: 0.5*np.cos(t)*(3/np.pi)**0.5
Y02 = lambda t: 0.25*((3*np.cos(t)*np.cos(t))-1)*(5/np.pi)**0.5

# get data
indices = []
xPos = []
yPos = []
zPos = []

# read data, populate lists
with open('disc.txt', mode='r') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter='\t')
    fields = next(csv_reader)
    for line in csv_reader:
        indices.append(int(line[0]))
        xPos.append(float(line[2]))
        yPos.append(float(line[3]))
        zPos.append(float(line[4]))
 
        
# get number of vertices for future use
Nv = len(indices)
print("Number of vertices: " + str(Nv))

# calculate derived values (spherical coordinates, radii) 
thetas = []
phis = []
radius = []
for i in range(0, Nv-1):
    radius.append(xPos[i]**2 + yPos[i]**2 + zPos[i]**2)
    if(zPos[i] == 0):
        phis.append(0)
    else:
        phis.append(np.arctan(((xPos[i])**2+(yPos[i])**2)**0.5/zPos[i]))
    if(xPos[i] == 0):
        thetas.append(0)
    else:
        thetas.append(np.arctan(yPos[i]/xPos[i]))
            
# create R vector, sampling from dataset (vertices 0, 100, 200)
R = []
avgR = sum(radius)/Nv
for i in range(0, 300, 100):
    R.append((radius[i] - avgR)/avgR)

# create square coefficient matrix, sampling from dataset (vertices 0, 100, 200)
S2 = np.zeros((3, 3))
row = 0
for i in range(0, 300, 100):
    S2[row, 0] = Y00
    S2[row, 1] = Y01(thetas[i])
    S2[row, 2] = Y02(thetas[i])
    row = row + 1

Alms = np.linalg.solve(S2, R)
print("Alms: ")
print(Alms)
