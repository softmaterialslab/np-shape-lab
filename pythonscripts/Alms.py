import numpy as np
import csv
import scipy.linalg as la

# define spherical harmonics constants
Y00 = 0.5*(1/np.pi)**0.5
Y01 = lambda t: 0.5*np.cos(t)*(3/np.pi)**0.5
Y02 = lambda t: 0.25*((3*np.cos(t)*np.cos(t))-1)*(5/np.pi)**0.5

# get data
indices = []
thetas = []
phis = []

with open('poles.txt', mode='r') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter='\t')
    fields = next(csv_reader)
    i=0
    for line in csv_reader:
        indices.append(int(line[0]))
        thetas.append(float(line[1]))
        phis.append(float(line[2]))
        i=i+1
        
# get number of vertices for future use
Nv = i
print("Number of vertices: " + str(Nv))
        
# calculate center of mass
COM = [0,0]
COM[0] = sum(thetas) / len(thetas)
COM[1] = sum(phis) / len(phis)

S = np.array([[Y00, Y01(thetas[0])],[Y00,Y01(thetas[1])]])

# TODO: calculate radius, r, so need to get cartesian coordinates
R = [[4],[4]]

x = np.linalg.solve(S, R)
print(x)


print(S)

