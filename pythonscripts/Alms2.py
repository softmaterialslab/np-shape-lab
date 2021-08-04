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
        
# calculate center of mass
COM = [0,0,0]
COM[0] = sum(xPos) / len(xPos)
COM[1] = sum(yPos) / len(yPos)
COM[2] = sum(zPos) / len(zPos)
print("Center of Mass:" + str(COM))

# calculate derived values (spherical coordinates, radii) 
thetas = []
phis = []
radius = []
for i in range(0, Nv-1):
    thetas.append(np.arctan(yPos[i]/xPos[i]))
    phis.append(np.arctan(((xPos[i])**2+(yPos[i])**2)**0.5/zPos[i]))
    radius.append(xPos[i]**2 + yPos[i]**2 + zPos[i]**2)

# create and populate numpy array
S = np.zeros((Nv, 3))
for i in range(0, Nv-1):
    S[i, 0] = Y00
    S[i, 1] = Y01(thetas[i])
    S[i, 2] = Y02(thetas[i])


print(S)
    
# create R vector
R = []
avgR = sum(indices)/Nv
for i in range(0, Nv-1):
    R.append((radius[i] - avgR)/avgR)


# use linalg to solve for Alms
#Alms = np.linalg.solve(S, R)
#print(Alms)


# testing: fill with zeros
T = np.zeros((Nv-1, Nv-1))
for i in range(0, Nv-1):
    T[i, 0] = Y00
    T[i, 1] = Y01(thetas[i])
    T[i, 2] = Y02(thetas[i])

Alms = np.linalg.solve(T,R)
