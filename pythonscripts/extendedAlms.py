import numpy as np
import csv
import scipy.linalg as la
import matplotlib.pyplot as plt

# indicate data file, should read from a p.lammpstrj file
filename = 'p.lammpstrj'

# define spherical harmonics constants
Y00 = 0.5*(1/np.pi)**0.5
Y10 = lambda t: 0.5*np.cos(t)*(3/np.pi)**0.5
Y20 = lambda t: 0.25*((3*np.cos(t)*np.cos(t))-1)*(5/np.pi)**0.5

# read data, populate lists
with open(filename, mode='r') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter='\t')

    # determine how many frames there are (total lines/381)
    frames = int(sum(1 for row in csv_reader) / 381)

    # reset file pointer to beginning
    csv_file.seek(0)

    # create arrays for data
    xPos = np.zeros((frames, 372))
    yPos = np.zeros((frames, 372))
    zPos = np.zeros((frames, 372))
    times = []

    # get data for each frame
    for frame in range(0,frames):

        # skip the first line, then get the timestep of the frame 
        next(csv_reader)
        times.append(next(csv_reader)[0])
        
        # skip the next 7 lines of each frame
        for i in range(0,7):
            next(csv_reader)

        # read in the next 371 lines (vertices)
        for j in range(0,372):
            line = next(csv_reader)
            xPos[frame,j]=line[2]
            yPos[frame,j]=line[3]
            zPos[frame,j]=line[4]

# for each frame, calculate derived values 
thetas = np.zeros((frames, 372))
phis = np.zeros((frames, 372))
radii = np.zeros((frames, 372))
averages = []

for frame in range(0,frames):
    
    for vertex in range(0, 372):
        radii[frame, vertex] = (xPos[frame, vertex]**2 + yPos[frame, vertex]**2 + zPos[frame, vertex]**2)
        if(zPos[frame, vertex] == 0):
            phis[frame, vertex] = 0
        else:
            phis[frame, vertex] = np.arctan(((xPos[frame, vertex])**2+(yPos[frame, vertex])**2)**0.5/zPos[frame, vertex])
        if(xPos[frame, vertex] == 0):
            thetas[frame, vertex] = 0
        else:
            thetas[frame, vertex] = np.arctan(yPos[frame, vertex]/xPos[frame, vertex])
    averages.append(sum(radii[frame])/372)
    
# calculate Alms based on sampled data (vertices 0, 124, 248)
Alms = np.zeros((frames, 3))
for frame in range(0, frames):
    # make S coefficient matrix and R vector
    S = np.zeros((3,3))
    R = []
    row = 0
    for i in range (0,372,124):
        S[row, 0] = Y00
        S[row, 1] = Y10(thetas[frame, i])
        S[row, 2] = Y20(thetas[frame, i])
        R.append((radii[frame, i] - averages[frame])/averages[frame])
        row = row + 1
    Alms[frame] = np.linalg.solve(S, R)
    print("Alms for Frame " + str(frame) + ": " + str(Alms[frame]))


# extract Alm data to plot
A0 = []
A1 = []
A2 = []
for i in range(0, frames):
    A0.append(Alms[i][0])
    A1.append(Alms[i][1])
    A2.append(Alms[i][2])
    
# graph Alms over time
plt.plot(times,A0,label='Y00')
plt.plot(times,A1,label='Y10')
plt.plot(times,A2,label='Y20')
plt.legend()
plt.title("Alm Values For " + filename + " (disc)")
plt.show()
