import Generation as gen
import numpy as np
import os
import sys
from ctypes import *
from QBladeLibrary import QBladeLibrary

#!!!!Cq contour broken since all data sent through transfer file and not transfer_v

def sim(v, x, airfoil, optimiser, iteration,): #x is always the input variables from the optimisation function. First is the chords, then the twsit then the offset. v is the car speed

    #split x variable into chord and twist and pos
    c = np.array([x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8],x[9]])
    theta = np.array([x[10],x[11],x[12],x[13],x[14],x[15],x[16],x[17],x[18],x[19]])

    netWind = 10 + v # v is car speed. 10 is average windspeed but this can be adjusted
    #Qblade library set up
    dll_directory = os.path.abspath("./QBlade")
    dll_files = [f for f in os.listdir(dll_directory) if 'QBlade' in f and ('.dll' in f or '.so' in f)]

    if not dll_files:
        print('No matching QBlade*.dll or QBlade*.so files found in the specified directory:',os.path.abspath(dll_directory))
        sys.exit(1)  # Exit the script with a non-zero status to indicate an error

    dll_file_path = os.path.join(dll_directory, dll_files[0])
    QBLIB = QBladeLibrary(dll_file_path)  

    QBLIB.createInstance(1,32)

    bladelength = 1

    cq = np.array([])
    ct = np.array([])
    tsr = np.array([])

    for i in range(0,25): # i / 2 is tsr
        w = i / 2 * netWind / bladelength
        rpm = w / 2 / np.pi * 60
        gen.createBlade(c, theta, airfoil, w, v)
        gen.createSim(netWind, rpm)

        QBLIB.loadSimDefinition(b"./Generated/sim.txt") 

        QBLIB.initializeSimulation()
        for y in range(150):
            QBLIB.advanceTurbineSimulation()

        QBLIB.exportResults(1,b"./sim_data",f"results_{i / 2}".encode(),b"")

        output = np.loadtxt(f"./sim_data/results_{i / 2}.dat")
        if np.isnan(output[16]):
            output[16] = 0

        cq = np.append(cq, [output[15]])

        ct = np.append(ct, [output[17]])

        tsr = np.append(tsr, [i / 2])

    combined =  np.column_stack((tsr, cq, ct))
    file_path = f"./transfers/{optimiser}/{iteration}.txt"
    os.makedirs(os.path.dirname(file_path), exist_ok= True)
    transfer = np.savetxt(file_path, combined, header="time, cq, ct", delimiter=" ")

    QBLIB.closeInstance()

    del QBLIB.lib

v = np.array([0,1,2,3,4,5])
x = np.array([0.1, 0.085, 0.07, 0.055, 0.04, 0.02, 0.02, 0.02, 0.01, 0, 0,15,10,10,10,5, 5, 0,0,0, 0,0,0.00,0.00,0.00])
#for i in v:
    #sim(i, x, "S1223")
# sim(4, [1.11383193e-01,5.67988936e-02,4.16924422e-02,5.45674496e-02,3.50018526e-02,3.80703120e-02,2.79462992e-02,2.73882467e-02,1.08476162e-02,1.02152241e-02,3.73000593e+01,2.89312871e+01,1.94715452e+01,1.98942402e+01,1.36604940e+01,1.37956521e+01,1.10388084e+01,8.81236874e+00,3.65453291e-01,5.48863049e-01,1.78609131e+00,5.09401003e+00,9.68859783e+00,2.29568822e+00,1.04270822e+00], "S1223", "aa6a76", "final")

sim(4, x, "S1223", "olszak", "0")