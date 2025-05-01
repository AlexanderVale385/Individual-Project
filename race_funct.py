import numpy as np
import scipy.integrate as scint
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import csv
import os

mass = 125 #kg
r_wheel = 0.1 #radius of the wheel
windspeed = -10 #Estimate. Quite high
Cr = 0.01 #assumed rolling resistance coefficient. From GDP
cD = 0.13 #From CFD in GDP
FrontalArea = 1.08 # From GDP
wheelefficiency = 0.95
gearefficiency = 0.9
rho = 1.225
rBlade = 1
rolling_resistance = Cr * mass * 9.81

#punishing variable global definition
deltawachrace = 0
punishingFactor = 0.01

#define gear change variable. Changed in 500m race file and changed in optimiser
vc1, vc2, g1, g2, g3 = (2.75,4,9,6,4)

#arrays to store data to be pushed to log file

#Reading output from Qblade
tsrQblade = cqQblade = ctQblade = cq_interp = ct_interp = None

#global variables for gear shifting
shifting = False
current_gear = 0
shift_duration = 1 #seconds
t_shift_start = None
next_gear = None
gear3Bool = False
def gear_ratio_func(t, v, vc1, vc2, gear_ratios):
    global shifting, current_gear, t_shift_start, next_gear

    # Check if shift should begin
    if not shifting:
        if current_gear ==  0 and v > vc1:
            next_gear = 1
            t_shift_start = t
            shifting = True
        elif current_gear == 1:
            if v > vc2:
                next_gear = 2
                t_shift_start = t
                shifting = True
            elif v < vc1 - 1:
                next_gear = 0
                t_shift_start = t
                shifting = True
        elif current_gear == 2 and v < vc2 - 1:
            next_gear = 1
            t_shift_start = t
            shifting = True

    # If shifting, blend gear ratios
    if shifting:
        dt = t - t_shift_start
        alpha = np.clip(dt / shift_duration, 0.0, 1.0)
        ratio = (1 - alpha) * gear_ratios[current_gear] + alpha * gear_ratios[next_gear]

        # Finish shift
        if alpha >= 1.0:
            current_gear = next_gear
            next_gear = None
            t_shift_start = None
            shifting = False

        return ratio
    return gear_ratios[current_gear]

# def ratio(v, vc1, vc2, g1, g2, g3):
#     if v > vc1 and v < vc2:
#         return g2
#     elif v > vc2:
#         return g3
#     else:
#         return g1
count = 0

def force(ratio, tsr, v):
    global deltawachrace, cq_interp, ct_interp

    vrel = v - windspeed

    minProp = (np.min(cqQblade) * 0.5 * rho * (vrel) ** 2 * np.pi * rBlade ** 3) * ratio / r_wheel * wheelefficiency * gearefficiency
    maxDrag = (np.max(ctQblade) * 0.5 * rho * (vrel) ** 2 * np.pi * rBlade ** 2) + rolling_resistance + (0.5 * rho * (vrel) ** 2 * cD * FrontalArea)
    if tsr > 12:
        #print(f"TSR high: {tsr}")
        deltawachrace += (tsr - 12) * punishingFactor #2 is punishing factor
        return minProp - maxDrag, minProp, maxDrag 
    elif tsr < 0:
        minDrag = (np.min(ctQblade) * 0.5 * rho * (vrel) ** 2 * np.pi * rBlade ** 2) + rolling_resistance + (0.5 * 1.225 * (vrel) ** 2 * cD * FrontalArea)
        #print(f"TSR low: {tsr}")
        deltawachrace -= tsr * punishingFactor #tsr is negative so 2 negatives make positive
        return minProp - minDrag, minProp, minDrag
    cq = cq_interp(tsr)
    ct = ct_interp(tsr)

    propulsion = (cq * 0.5 * rho * (vrel) ** 2 * np.pi * rBlade ** 3) * ratio / r_wheel * gearefficiency * wheelefficiency #np.pi is area of blade as r = 1.

    drag = (ct * 0.5 * rho * (vrel) ** 2 * np.pi * rBlade ** 2)  + rolling_resistance + (0.5 * rho * (vrel) ** 2 * cD * FrontalArea)

    net = propulsion - drag

    return net, propulsion, drag

def tsrCalc(ratio, v):
    omega_wheel = v / r_wheel
    tip_vel = omega_wheel * ratio  * rBlade

    tsr = tip_vel / (v - windspeed)

    return tsr

def acceleration(ratio, v): #v represents car speed
    tsr = tsrCalc(ratio, v)
    
    net, propulsion, drag = force(ratio, tsr, v)
    return net / mass, propulsion, drag

def ratio_velocity_force(optimiser, iteration):
    #loading correct data into program
    global tsrQblade, cqQblade, ctQblade, cq_interp, ct_interp
    output = np.loadtxt(f"./transfers/{optimiser}/{iteration}.txt", skiprows= 1)
    #Reading output from Qblade
    tsrQblade = output[:,0]
    cqQblade = output[:,1]
    ctQblade = output[:,2]

    cq_interp = interp1d(tsrQblade, cqQblade, kind = "linear")
    ct_interp = interp1d(tsrQblade, ctQblade, kind = "linear")

    max_line = np.empty((0, 3))
    data = np.empty((151,151))
    tsr = np.empty((151,151))
    x = np.array([])
    for y in range(151): #iterate through vel
        x = np.append(x, y / 10)
        yg = np.array([])
        z = np.array([])

        for i in range(151): #iterate through ratios
            a, p, d = acceleration(i / 10, y / 10)
            net = a * mass
            yg = np.append(yg, i / 10)
            z = np.append(z, net)
            data[y, i] = net
            tsr[y, i] = tsrCalc(i / 10, y / 10)
            
        net_max = np.max(z)
        rat_max = yg[np.where(z == net_max)[0][0]]
        max_line = np.vstack([max_line, np.array([y / 10, rat_max, net_max])])

    return data, max_line, x, yg, tsr

def ratio_velocity_prop(optimiser, iteration):
    #loading correct data into program
    global tsrQblade, cqQblade, ctQblade, cq_interp, ct_interp
    output = np.loadtxt(f"./transfers/{optimiser}/{iteration}.txt", skiprows= 1)
    #Reading output from Qblade
    tsrQblade = output[:,0]
    cqQblade = output[:,1]
    ctQblade = output[:,2]

    cq_interp = interp1d(tsrQblade, cqQblade, kind = "linear")
    ct_interp = interp1d(tsrQblade, ctQblade, kind = "linear")


    max_line = np.empty((0, 3))
    data = np.empty((81,81))
    x = np.array([])
    tsr = np.empty((81,81))
    for i in range(81):
        x = np.append(x, i / 10)
        yv = np.array([])
        z = np.array([])

        for y in range(81):
            f, prop, d = acceleration(i / 10, y / 10)
            yv = np.append(yv, y / 10)
            z = np.append(z, prop)
            data[y, i] = prop
            tsr[y, i] = tsrCalc(i / 10, y / 10)

        prop_max = np.max(z)
        vel_max = yv[np.where(z == prop_max)[0][0]]
        max_line = np.vstack([max_line, np.array([i / 10, vel_max, prop_max])])

    return data, max_line, x, yv, tsr

def ratio_velocity_drag(optimiser, iteration):
    #loading correct data into program
    global tsrQblade, cqQblade, ctQblade, cq_interp, ct_interp
    output = np.loadtxt(f"./transfers/{optimiser}/{iteration}.txt", skiprows= 1)
    #Reading output from Qblade
    tsrQblade = output[:,0]
    cqQblade = output[:,1]
    ctQblade = output[:,2]

    cq_interp = interp1d(tsrQblade, cqQblade, kind = "linear")
    ct_interp = interp1d(tsrQblade, ctQblade, kind = "linear")

    max_line = np.empty((0, 3))
    data = np.empty((151,151))
    x = np.array([])
    tsr = np.empty((151,151))
    for i in range(151):
        x = np.append(x, i / 10)
        yv = np.array([])
        z = np.array([])

        for y in range(151):
            f, prop, drag = acceleration(i / 10, y / 10)
            yv = np.append(yv, y / 10)
            z = np.append(z, drag)
            data[y, i] = drag
            tsr[y, i] = tsrCalc(i / 10, y / 10)

        drag_max = np.max(z)
        vel_max = yv[np.where(z == drag_max)[0][0]]
        max_line = np.vstack([max_line, np.array([i / 10, vel_max, drag_max])])

    return data, max_line, x, yv,tsr

def model(t, f, vc1, vc2, g1, g2, g3):
    global deltawachrace, gear3Bool
    x = f[0]
    v = f[1]
    gear_ratio = gear_ratio_func(t, v, vc1, vc2, [g1, g2, g3])
    if gear_ratio == g3:
        gear3Bool = True
    a, prop, drag = acceleration(gear_ratio, v)
    tsr = tsrCalc(gear_ratio, v)
    #print(f"t={t:.2f}, v={v:.2f}, a={a:.4f}, netF={prop - drag:.2f}, prop={prop:.2f}, drag={drag:.2f}, tsr = {tsr}, x = {x}")
    return [v, a]

delta_t = 0.01

def race(optimiser, iteration, v0, vc1, vc2, g1, g2, g3):
    global deltawachrace, tsrQblade, cqQblade, ctQblade, cq_interp, ct_interp, gear3Bool
    gear3Bool = False
    y0 = (0, v0)
    deltawachrace = 0

    output = np.loadtxt(f"./transfers/{optimiser}/{iteration}.txt", skiprows= 1)
    #Reading output from Qblade
    tsrQblade = output[:,0]
    cqQblade = output[:,1]
    ctQblade = output[:,2]

    cq_interp = interp1d(tsrQblade, cqQblade, kind = "linear")
    ct_interp = interp1d(tsrQblade, ctQblade, kind = "linear")

    def finish_race(t, y, vc1, vc2, g1, g2, g3):
        return y[0] - 600
    
    def finish_accel(t, y, vc1, vc2, g1, g2, g3):
        return y[0] - 100
    
    finish_race.terminal = True
    finish_accel.direction = 1
    finish_race.direction = 1

    out = scint.solve_ivp(model, (0,1000), y0, events = [finish_accel, finish_race], dense_output= True, args = (vc1, vc2, g1, g2, g3), method = "BDF", max_step = delta_t)
    
    if out.success == False:
        return [-10000000], None, out, deltawachrace
    
    t_finish_accel= out.t_events[0]
    t_finish_race = out.t_events[1]

    average_speed = 500 / (t_finish_race - t_finish_accel)

    wach = average_speed / -windspeed

    vend = out.y[1][np.where(out.t == t_finish_race)]

    print(f"Gear 3: {gear3Bool}")
    if gear3Bool == False:
        deltawachrace += (100 / vend)

    return wach, vend, out, deltawachrace

def log(optimiser, iteration, out, x, extra, wach, error, bladeOpt = bool):
    global cq_interp, ct_interp
    tsrArray = np.array([])
    cqArray = np.array([])
    ctArray = np.array([])
    dragArray = np.array([])
    propArray = np.array([])
    netArray = np.array([])
    directory = f"./logs/{optimiser}"
    os.makedirs(directory, exist_ok=True)
    # Define your columns (headers)
    columns = ['X', 'Wach', 'Error', 'Iteration']

    # File path

    if bladeOpt == True:
        destination = f"./logs/{optimiser}/mainBlade.csv"
    else:
        destination = f"./logs/{optimiser}/mainGear.csv"
    # Check if the file already exists
    file_exists = os.path.exists(destination)
    with open(destination, "a", newline = '') as file:
        log = csv.writer(file)

        # Write header only if the file is empty
        if not file_exists:
            log.writerow(columns)  # Write the header row

        totalTime = out.t
        totalVel = out.y[1]
        totalDis = out.y[0]
        for i in range(len(totalTime)):
            #change which objective variables to use depending on if it is blade log or gear log
            if bladeOpt == True:
                gear = gear_ratio_func(totalTime[i], totalVel[i], extra[0], extra[1], [extra[2], extra[3], x[24]])
            else:
                gear = gear_ratio_func(totalTime[i], totalVel[i], x[0], x[1], [x[2], x[3], extra[0]])
            tsr = tsrCalc(gear, totalVel[i])

            try:
                cq = cq_interp(tsr)
                ct = ct_interp(tsr)
            except ValueError:
                cq = 0
                ct = 0

            propulsion = (cq * 0.5 * rho * (totalVel[i]-windspeed) ** 2 * np.pi * rBlade ** 3) * gear / r_wheel * gearefficiency * wheelefficiency#no-slip assumed. np.pi is area of blade as r = 1.

            drag = (ct * 0.5 * rho * (totalVel[i] - windspeed) ** 2 * np.pi * rBlade ** 2)  + rolling_resistance + (0.5 * rho * (totalVel[i] - windspeed) ** 2 * cD * FrontalArea)

            net = propulsion - drag

            tsrArray = np.append(tsrArray, tsr)
            cqArray = np.append(cqArray, cq)
            ctArray = np.append(ctArray, ct)
            dragArray = np.append(dragArray, drag)
            propArray = np.append(propArray, propulsion)
            netArray = np.append(netArray, net)

        arrays = [totalTime, totalVel, tsrArray, cqArray, ctArray, dragArray, propArray, netArray, totalDis]
        np.savetxt(f"./logs/{optimiser}/{iteration}.csv", arrays, delimiter = ",", fmt = "%f")



        log.writerows([[x, wach, error, iteration]])
iterationD = 0
def standardRaceGearOpt(x, optimiser, referenceOpt, referenceIt):
    global iterationD
    deltawach = 0
    #fix in order to stop 2 gear change solution and failure upon increasing gear ratio at high speeds.
    if x[3] < x[4]:
        deltawach = (x[4] - x[3]) * 100
    if x[2] < x[3]:
        deltawach = (x[3] - x[2]) * 100
    x = np.array(x)
    print(x)
    #5m push start gives initial velocity of 1 m/s. For test blade, this is 3m/s to make it work
    wach, vend, out, delta1 = race(referenceOpt, referenceIt, 1, x[0], x[1], x[2], x[3], x[4]) # acceleration period of race. v1 is speed at start of race time. w is not used

    print(wach, vend)
    deltawach += delta1
    #comment out to disable logging function
    log(optimiser, iterationD + referenceIt, out, x, [x[4]], wach, deltawach, False)

    #Punish for straying
    print(deltawach)
    wach = wach - deltawach

    iterationD += 1
    return 0 - wach



