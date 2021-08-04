from numpy import linalg as LA
import csv
import matplotlib.pyplot as plt
import numpy as np

indices = []
xPos = []
yPos = []
zPos = []
charge = []

# loop through data as a csv, store data in lists
with open('sphere.txt', mode='r') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter='\t')
    fields = next(csv_reader)
    i=0
    for line in csv_reader:
        indices.append(int(line[0]))
        xPos.append(float(line[2]))
        yPos.append(float(line[3]))
        zPos.append(float(line[4]))
        charge.append(float(line[5]))
        i=i+1
      
# get number of vertices for future use
Nv = i
print("Number of vertices: " + str(Nv))

# calculate center of mass
COM = [0,0,0]
COM[0] = sum(xPos) / len(xPos)
COM[1] = sum(yPos) / len(yPos)
COM[2] = sum(zPos) / len(zPos)
print("Center of Mass:" + str(COM))

# set up matrix S
S = [[0 for x in range(3)] for y in range(3)]
diag1, diag2, diag3, xy, xz, yz = 0, 0, 0, 0, 0, 0
# populate diagonals and off-diagonals
for x in range(0, Nv):
    diag1 = diag1 + (xPos[x]-COM[0])**2
    diag2 = diag2 + (yPos[x]-COM[1])**2
    diag3 = diag3 + (zPos[x]-COM[2])**2

    xy = xy + (xPos[x]-COM[0])*(yPos[x]-COM[1])
    xz = xz + (xPos[x]-COM[0])*(zPos[x]-COM[2])
    yz = yz + (yPos[x]-COM[1])*(zPos[x]-COM[2])
    
S[0][0] = diag1 / Nv
S[0][1] = xy / Nv
S[0][2] = xz / Nv
S[1][0] = xy / Nv
S[1][1] = diag2 / Nv
S[1][2] = yz / Nv
S[2][0] = xz / Nv
S[2][1] = yz / Nv
S[2][2] = diag3 / Nv

# calculate eigenvalues
S = np.array(S)
w, v = LA.eig(S)
print("Eigenvalues of S: " + str(w))

# calculate squared radius of gyration
Rg2 = (w[0])**2 + (w[1])**2 + (w[2])**2

# calculate normalized asphericitty (equation 2 of 2020 paper)
asphericity = ((1.5*((max(w))**2)) - (0.5*(Rg2))) / Rg2

print("Asphericity: " + str(asphericity))

# graph deviations

# calculate magnitude of distances
distance = []
for i in range(0, Nv):
    distance.append((xPos[i]**2)+(yPos[i]**2)+(zPos[i]))
    if distance[i] < 0:
        distance[i] = -1 * distance[i]
indices = np.arange(0,Nv)
avgRadius = []
sphere = np.arange(0,Nv)
for item in sphere:
    sphere[item]=1

average1 = sum(distance)/Nv

for item in range(0,Nv):
    avgRadius.append(average1)


plt.plot(indices, distance, 'o', label='vertices distance')
plt.plot(indices, avgRadius, label='average distance')
plt.plot(indices,sphere,label='original distance')
plt.legend()
plt.title("Disc Dataset")
plt.show()
