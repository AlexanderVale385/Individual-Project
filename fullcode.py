import QbladeWrapper as Qblade
import race_funct as race
import scipy.optimize as sco
import numpy as np
from ctypes import *
import uuid
import csv
import pandas as pd


solve = uuid.uuid4().hex[:6]
solve = "a6aa76"
iterationB = 1501

def standardRaceBldOpt(x):
    global deltawach, iterationB
    deltawach = 0
    Qblade.sim(4, x, airfoil, solve, iterationB)

    # minChord = 0.01 #chord cant go below this value
    # maxChord = 0.15
    if x[23] < x[24]:
        print("gear change 2 called")
        deltawach += (x[24] - x[23]) * 100
    if x[22] < x[23]:
        print("gear change 1 called")
        deltawach += (x[23] - x[22]) * 100

    
    #fix in order to stop 2 gear change solution and failure upon increasing gear ratio at high speeds.
    x = np.array(x)
    print(x)
    #5m push start gives initial velocity of 1 m/s. For test blade, this is 3m/s to make it work
    wach, vend, out, delta1 = race.race(solve, iterationB, 1, x[20], x[21], x[22], x[23], x[24]) # acceleration period of race. v1 is speed at start of race time. w is not used

    print(wach, vend)
    deltawach = deltawach + delta1
    #comment out to disable logging function
    race.log(solve, iterationB, out, x, [x[20], x[21], x[22], x[23]], wach, deltawach, True)

    #Punish for straying
    print(deltawach)
    wach = wach - deltawach
    iterationB += 1
    return 0 - wach

#inital simplex set up:
airfoil = "S1223"
g0 = [2.1,3.9,6.2,3.4]
g30 = 1.4
x0 = np.array([0.1, 0.074, 0.053, 0.042, 0.034, 0.029, 0.025, 0.022, 0.02, 0.018, 32.5, 18.6, 12.1, 8.4, 6, 4.4, 3.2, 2.2, 1.5, 0.9, g0[0], g0[1], g0[2], g0[3], g30]) #last 1 parameters are gear optimisation
initial_simplex = [x0]
ChordBounds, TwistBounds = [[(0.01, 0.15)], [(0, 80)]]
geometry_param_total = (len(x0) - 5) / 2
Bladebounds = []
chord_simplex_grad = -0.1 / geometry_param_total
theta_simplex_grad = -5 / geometry_param_total

for i in range(len(x0)):
    y = x0.copy()
    if i < geometry_param_total:
        Bladebounds.extend(ChordBounds)
        z = 0.01 #inital simplex for chord
    elif i < 2 * geometry_param_total:
        Bladebounds.extend(TwistBounds)
        z = 10 #initial simplex for twist
    elif i < 2 * geometry_param_total + 2:
        Bladebounds.extend([(1, 12)])
        z = -1
    else:
        Bladebounds.extend([(0.1, 10)])
        z = -1

    y[i] += z

    initial_simplex.append(y)
print(Bladebounds)

def callbackB(xk):
    print(f"Callback called: Current best point: {xk}")
    with open(f"./logs/control points/{solve}_Blade.csv", "a", newline = '') as f:
        writer = csv.writer(f)
        writer.writerow(xk.tolist())
def callbackG(xk):
    print(f"Callback called: Current best point: {xk}")
    with open(f"./logs/control points/{solve}_Gear.csv", "a", newline = '') as f:
        writer = csv.writer(f)
        writer.writerow(xk.tolist())

def manage(word):
    if type(word) == float or type(word) == int:
        return word
    array = np.array(word.split(" "))
    remove_chars = "array()/n[]"
    clean_func = np.vectorize(lambda s: s.translate(str.maketrans("", "", remove_chars)))
    clean = clean_func(array)
    cleaner = clean[clean != ""]
    cleanest = cleaner[cleaner != "..."]
    arrayint = cleanest.astype(float)
    return np.array(arrayint)

# print(standardRaceBldOpt([1.11383193e-01,5.67988936e-02,4.16924422e-02,5.45674496e-02,3.50018526e-02,3.80703120e-02,2.79462992e-02,2.73882467e-02,1.08476162e-02,1.02152241e-02,3.73000593e+01,2.89312871e+01,1.94715452e+01,1.98942402e+01,1.36604940e+01,1.37956521e+01,1.10388084e+01,8.81236874e+00,3.65453291e-01,5.48863049e-01,1.78609131e+00,5.09401003e+00,9.68859783e+00,2.29568822e+00,1.04270822e+00]))
Qblade.sim(4, [1.11383193e-01,5.67988936e-02,4.16924422e-02,5.45674496e-02,3.50018526e-02,3.80703120e-02,2.79462992e-02,2.73882467e-02,1.08476162e-02,1.02152241e-02,3.73000593e+01,2.89312871e+01,1.94715452e+01,1.98942402e+01,1.36604940e+01,1.37956521e+01,1.10388084e+01,8.81236874e+00,3.65453291e-01,5.48863049e-01,1.78609131e+00,5.09401003e+00,9.68859783e+00,2.29568822e+00,1.04270822e+00], airfoil, "aa6a76", "final")

# solve = "12b6ac"
# AccelResult1 = sco.minimize(race.standardRaceGearOpt, g0, args = (g30, solve, solve, iterationB),options = {"maxfev": 100}, method = "Nelder-Mead", tol = 1e-4, callback=callbackG)
# main = pd.read_csv(f'./logs/{solve}/mainGear.csv', skiprows=0).to_numpy()
# data = np.array([[manage(word) for word in row] for row in main]) #0: x, 1: wach, 2:error, 3: iterationB

#looking for this condition
# maxWach = np.max(data[:, 1])
# for i in range(len(data[:, 0])):
#     if data[i, 1] == maxWach:
#         vc1, vc2, g1, g2 = data[i, 0]

# iterationB = np.max(data[:, 3]) + 1

# res = sco.minimize(standardRaceBldOpt, x0, method = "Nelder-Mead", tol = 1e-4, options = {"initial_simplex": initial_simplex}, callback= callbackB, bounds= Bladebounds)
# print(res) 

# main2 = pd.read_csv(f'./logs/{solve}/mainBlade.csv', skiprows=0).to_numpy()
# data2 = np.array([[manage(word) for word in row] for row in main2]) #0: x, 1: wach, 2:error, 3: iterationB

#looking for this condition
# maxWach2 = np.max(data2[:, 1])
# for i in range(len(data2[:, 0])):
#     if data2[i, 1] == maxWach2:
#         x1 = data2[i, 0]
#         iterationData = data2[i, 3]

#AccelResult2 = sco.minimize(race.standardRaceGearOpt, [vc1, vc2, g1, g2], args = (x1[20], solve, solve, iterationData), method = "Nelder-Mead", tol = 1e-4, callback=callbackG)

