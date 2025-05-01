#Set up template scripts
import numpy as np
import bisect
from scipy.interpolate import make_interp_spline
import time

trb = """----------------------------------------QBlade Turbine Definition File----------------------------------------------
Generated with : QBlade CE v2.0.7.7_beta windows
Archive Format: 310024
Time : 17:36:33
Date : 14.10.2024

----------------------------------------Object Name-----------------------------------------------------------------
generated                             OBJECTNAME         - the name of the turbine object

----------------------------------------Rotor Definition------------------------------------------------------------
bld.txt                        BLADEFILE          - the path of the blade file that is used in this turbine definition
0                                        TURBTYPE           - the turbine type (0 = HAWT or 1 = VAWT)
3                                        NUMBLADES          - the number of blades (a Structural Model overrides this value)
0                                        ROTORCONFIG        - the rotor configuration (0 = UPWIND or 1 = DOWNWIND)
0                                        ROTATIONALDIR      - the direction of rotor rotation (0 = STANDARD or 1 = REVERSED)
2                                        DISCTYPE           - type of rotor discretization (0 = from Bladetable, 1 = linear, 2 = cosine) 
0                                        INTPTYPE           - type of rotor interpolation (0 = linear, 1 = Akima splines) 
150                                       NUMPANELS          - the number of aerodynamic panels per blade (unused if DISCTYPE = 0)

----------------------------------------Turbine Geometry Parameters-------------------------------------------------
These values are only used if no Structural Model is defined for this Turbine, in case of a Structural Model the geometry is defined in the Structural Input Files!!
0.1670                                  OVERHANG           - the rotor overhang [m] (HAWT only)
0.0000                                   SHAFTTILT          - the shaft tilt angle [deg] (HAWT only)
0.0000                                   ROTORCONE          - the rotor cone angle [deg] (HAWT only)
0.0000                                   CLEARANCE          - the rotor clearance to ground [m] (VAWT only)
0.0000                                   XTILT              - the rotor x-tilt angle [deg] (VAWT only)
0.0000                                   YTILT              - the rotor y-tilt angle [deg] (VAWT only)
2.0000                                 TOWERHEIGHT        - the tower heigh [m]
0.029                                  TOWERTOPRAD        - the tower top radius [m]
0.04                                  TOWERBOTRAD        - the tower bottom radius [m]

----------------------------------------Dynamic Stall Models--------------------------------------------------------
0                                        DYNSTALLTYPE       - the dynamic stall model: 0 = none; 1 = OYE; ; 2 = IAG; 3 = GORMONT-BERG; 4 = ATEFLAP
8                                        TF_OYE             - Tf constant for the OYE dynamic stall model
6                                        AM_GB              - Am constant for the GORMONT-BERG dynamic stall model
3                                        TF_ATE             - Tf constant for the ATEFLAP dynamic stall model
2                                        TP_ATE             - Tp constant for the ATEFLAP dynamic stall model
IAGPARAMS           - the parameters for the IAG DS model, order according to docs & GUI
0.30 0.70 0.70 0.53 0.75 1.70 3.00 6.00 6.00 0.20 0.10 1.50 1.50 0.00 

----------------------------------------Aerodynamic Models----------------------------------------------------------
true                                    UNSTEADYAERO       - include unsteady attached flow aerodynamics? [bool]
true                                     2PLIFTDRAG         - include the 2 point lift drag correction? [bool]
false                                    HIMMELSKAMP        - include the Himmelskamp Stall delay? (HAWT only) [bool]
false                                    TOWERSHADOW        - include the tower shadow effect [bool]
0.0                                     TOWERDRAG          - the tower drag coefficient [-] (if a Structural Model is used the tower drag is defined in the tower input file)

----------------------------------------Wake Type------------------------------------------------------------------
1                                        WAKETYPE           - the wake type: 0 = free vortex wake; 1 = unsteady BEM (unsteady BEM is only available for HAWT)

----------------------------------------Vortex Wake Parameters------------------------------------------------------
Only used if waketype = 0
0                                        WAKEINTTYPE        - the wake integration type: 0 = EF; 1 = PC; 2 = PC2B
true                                     WAKEROLLUP         - calculate wake self-induction [bool]
true                                     TRAILINGVORT       - include trailing vortex elements [bool]
true                                     SHEDVORT           - include shed vortex elements [bool]
0                                        CONVECTIONTYPE     - the wake convection type (0 = BL, 1 = HH, 2 = LOC)
1.00                                     WAKERELAXATION     - the wake relaxation factor [0-1]
1.00                                     FIRSTWAKEROW       - first wake row length [-]
200000                                   MAXWAKESIZE        - the maximum number of wake elements [-]
100                                      MAXWAKEDIST        - the maxmimum wake distance from the rotor plane (normalized by dia) [-]
0.05                                 WAKEREDUCTION      - the wake reduction factor [-]
0                                        WAKELENGTHTYPE     - the wake length type (0 = counted in rotor revolutions, 1 = counted in time steps)
1000000.00                               CONVERSIONLENGTH   - the wake conversion length (to particles) [-]
0.50                                     NEARWAKELENGTH     - the near wake length [-]
2.00                                     ZONE1LENGTH        - the wake zone 1 length [-]
4.00                                     ZONE2LENGTH        - the wake zone 2 length [-]
6.00                                     ZONE3LENGTH        - the wake zone 3 length [-]
2                                        ZONE1FACTOR        - the wake zone 1 factor (integer!) [-]
2                                        ZONE2FACTOR        - the wake zone 2 factor (integer!) [-]
2                                        ZONE3FACTOR        - the wake zone 3 factor (integer!) [-]

----------------------------------------Vortex Core Parameters------------------------------------------------------
Only used if waketype = 0
0.05                                     BOUNDCORERADIUS    - the fixed core radius of the bound blade vortex (fraction of local chord) [0-1]
0.05                                     WAKECORERADIUS     - the intial core radius of the free wake vortex (fraction of local chord) [0-1]
20                                   VORTEXVISCOSITY    - the turbulent vortex viscosity
false                                    VORTEXSTRAIN       - calculate vortex strain [bool]
50                                       MAXSTRAIN          - the maximum element strain, before elements are removed from the wake [-]

----------------------------------------Gamma Iteration Parameters--------------------------------------------------
Only used if waketype = 0
0.100                                    GAMMARELAXATION    - the relaxation factor used in the gamma (circulation) iteration [0-1]
0.001                                 GAMMAEPSILON       - the relative gamma (circulation) convergence criteria
100                                     GAMMAITERATIONS    - the maximum number of gamma (circulation) iterations (integer!) [-]

----------------------------------------Unsteady BEM Parameters------------------------------------------------------
12                                       POLARDISC          - the polar discretization for the unsteady BEM (integer!) [-]
true                                     BEMTIPLOSS         - use BEM tip loss factor [bool]
100                                      BEMSPEEDUP         - initial BEM convergence acceleration time [s]

----------------------------------------Structural Model-------------------------------------------------------------
			                 STRUCTURALFILE     - the input file for the structural model (leave blank if unused)
false                                    GEOMSTIFFNESS      - enable geometric stiffness [bool]
false                                    AEROPANELLOADS     - enable distributed aero panel loads and gradients [bool]

----------------------------------------Turbine Controller-----------------------------------------------------------
0                                        CONTROLLERTYPE     - the type of turbine controller 0 = none, 1 = BLADED, 2 = DTU, 3 = TUB
		                         CONTROLLERFILE     - the controller file name, WITHOUT file ending (.dll or .so ) - leave blank if unused
					   PARAMETERFILE      - the controller parameter file name (leave blank if unused)

TURB_INFO
v1.2-------DM------21.05.2024
now using Timoshenko instead of Euler-Bernoulli beam


END_TURB_INFO"""

def createBlade(c, theta, airfoil, w, v):
    pos = np.linspace(0.1, 1, len(c))
    x_fine = np.linspace(0.1, 1, 90)

    chordSpline = make_interp_spline(pos, c, k = 3)
    chord_fine = chordSpline(x_fine)

    thetaSpline = make_interp_spline(pos, theta, k = 3)
    theta_fine = thetaSpline(x_fine)

    if airfoil == "S1223":
        polarReynolds = [50, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000] #these are available reynolds numbers for this particular airfoil. In thousands

    reynolds = np.sqrt((x_fine * w) ** 2 + v ** 2) * 1.225 * chord_fine / (1.7894 * 10 ** -5) / 1000 #divide by 1000 to make it more legible in 1000s. Pythagoras applied to include calculation of oncoming air
    closestReynolds = np.array([], dtype= object)
    for i in reynolds:
        loc = bisect.bisect_left(polarReynolds, i)

        if loc == 0:
            closestReynolds = np.append(closestReynolds, polarReynolds[0])
        elif loc == len(polarReynolds):
            closestReynolds = np.append(closestReynolds, polarReynolds[-1])
        else:
            before = polarReynolds[loc - 1]
            after = polarReynolds[loc]
            if abs(before - i) <= abs(after - i):
                closestReynolds = np.append(closestReynolds, before)
            else:
                closestReynolds = np.append(closestReynolds, after)

    for i in range(90):
        if closestReynolds[i] == 50:
            closestReynolds[i] = "0.050"
        elif closestReynolds[i] == 1000:
            closestReynolds[i] = "1.000"
        else:
            closestReynolds[i] = f"0.{closestReynolds[i]}"
    
    airfoilarr = np.array([])

    for i in range(90):
        airfoilarr = np.append(airfoilarr, f"Polars/{airfoil}_Re{closestReynolds[i]}_M0.00_N9.0_360_M.plr")
    
    bld = f"""----------------------------------------QBlade Blade Definition File------------------------------------------------
Generated with : QBlade CE v2.0.7.7_beta windows
Archive Format: 310024
Time : 17:36:33
Date : 14.10.2024

----------------------------------------Object Name-----------------------------------------------------------------
generated                                 OBJECTNAME         - the name of the blade object

----------------------------------------Parameters------------------------------------------------------------------
HAWT                                     ROTORTYPE          - the rotor type
false                                    INVERTEDFOILS      - invert the airfoils? (only VAWT) [bool]
3                                        NUMBLADES          - number of blades

----------------------------------------Blade Data------------------------------------------------------------------
POS_[m]             CHORD_[m]           TWIST_[deg]         OFFSET_X_[m]        OFFSET_Y_[m]        P_AXIS [-]          POLAR_FILE  
0                    0.03                   0            0.00000             0.00000             0.25000               Polars/Circular_Foil_Circular_Foil_CD1.20.plr
0.05                    0.03                   0            0.00000             0.00000             0.25000               Polars/Circular_Foil_Circular_Foil_CD1.20.plr        
"""
    
    try:
        with open("./Generated/bld.txt", "w") as f:
                f.writelines(bld)
    except PermissionError:
        time.sleep(0.1)
        with open("./Generated/bld.txt", "w") as f:
                f.writelines(bld)
    

    try:
        with open("./Generated/bld.txt", "a") as f:
            for i in range(90):
                f.writelines(f"{x_fine[i]:.8f}                 {chord_fine[i]:.8f}             {theta_fine[i]:.8f}           0.00000             0.00000             0.25000            {airfoilarr[i]}\n")
    except PermissionError:
        time.sleep(0.1)
        with open("./Generated/bld.txt", "a") as f:
            for i in range(90):
                f.writelines(f"{x_fine[i]:.8f}                 {chord_fine[i]:.8f}             {theta_fine[i]:.8f}           0.00000             0.00000             0.25000            {airfoilarr[i]}\n")
    

def createSim(windSpeed, rpm):
    sim = f"""----------------------------------------QBlade Simulation Definition File------------------------------------------
Generated with : QBlade CE v2.0.7-release_candidate_beta windows
Archive Format: 310023
Time : 21:27:48
Date : 15.05.2024

----------------------------------------Object Name-----------------------------------------------------------------
generated                            OBJECTNAME         - the name of the simulation object

----------------------------------------Simulation Type-------------------------------------------------------------
0                                        ISOFFSHORE         - use a number: 0 = onshore; 1 = offshore

----------------------------------------Turbine Parameters---------------------------------------------------------
multiple turbines can be added by adding multiple definitions encapsulated with TURB_X and END_TURB_X, where X must start at 1

TURB_1
    trb.txt                     TURBFILE           - the turbine definition file that is used for this simulation
    NREL_5MW_OC4_SEMI_RWT                TURBNAME           - the (unique) name of the turbine in the simulation (results will appear under this name)
    0                                INITIAL_YAW        - the initial turbine yaw in [deg]
    0                                INITIAL_PITCH      - the initial collective blade pitch in [deg]
    0                                INITIAL_AZIMUTH    - the initial azimuthal rotor angle in [deg]
    1                                    STRSUBSTEP         - the number of structural substeps per timestep (usually 1)
    1                                    RELAXSTEPS         - the number of initial static structural relaxation steps
    1                                    PRESCRIBETYPE      - rotor RPM prescribe type (0 = ramp-up; 1 = whole sim; 2 = no RPM prescibed)
    {rpm}                                RPMPRESCRIBED      - the prescribed rotor RPM [-]
    5                                    STRITERATIONS      - number of iterations for the time integration (used when integrator is HHT or Euler)
    0                                    MODNEWTONITER      - use the modified newton iteration?
    1                                    INCLUDEAERO        - include aerodynamic forces?
    1                                    INCLUDEHYDRO       - include hydrodynamic forces?
    0.00                                 GLOBPOS_X          - the global x-position of the turbine [m]
    0.00                                 GLOBPOS_Y          - the global y-position of the turbine [m]
    0.00                                 GLOBPOS_Z          - the global z-position of the turbine [m]
    0.00                                 GLOBROT_X          - the global x-rotation of the turbine [deg]
    0.00                                 GLOBROT_Y          - the global y-rotation of the turbine [deg]
    0.00                                 GLOBROT_Z          - the global z-rotation of the turbine [deg]
                                         EVENTFILE          - the file containing fault event definitions (leave blank if unused)
                                         LOADINGFILE        - the loading file name (leave blank if unused)
                                         SIMFILE            - the simulation file name (leave blank if unused)
                                         MOTIONFILE         - the prescribed motion file name (leave blank if unused)
END_TURB_1

----------------------------------------Simulation Settings-------------------------------------------------------
0.00050000                                 TIMESTEP           - the timestep size in [s]
200                                        NUMTIMESTEPS       - the number of timesteps
0.000                                    RAMPUP             - the rampup time for the structural model
0.000                                    ADDDAMP            - the initial time with additional damping
50.000                                  ADDDAMPFACTOR      - for the additional damping time this factor is used to increase the damping of all components
0.000                                    WAKEINTERACTION    - in case of multi-turbine simulation the wake interaction start at? [s]

----------------------------------------Wind Input-----------------------------------------------------------------
0                                        WNDTYPE            - use a number: 0 = steady; 1 = windfield; 2 = hubheight
                                         WNDNAME            - filename of the turbsim input file, mann input file or hubheight file (with extension), leave blank if unused
0                                        STITCHINGTYPE      - the windfield stitching type; 0 = periodic; 1 = mirror
true                                     WINDAUTOSHIFT      - the windfield shifting automatically based on rotor diameter [bool]
0.00                                     SHIFTTIME          - the windfield is shifted by this time if WINDAUTOSHIFT = 0
{windSpeed}                                    MEANINF            - the mean inflow velocity, overridden if a windfield or hubheight file is use
0.00                                     HORANGLE           - the horizontal inflow angle
0.00                                     VERTANGLE          - the vertical inflow angle
0                                        PROFILETYPE        - the type of wind profile used (0 = Power Law; 1 = Logarithmic)
0.00                                    SHEAREXP           - the shear exponent if using a power law profile, if a windfield is used these values are used to calculate the mean wake convection velocities
0.010                                    ROUGHLENGTH        - the roughness length if using a log profile, if a windfield is used these values are used to calculate the mean wake convection velocities
0.00                                     DIRSHEAR           - a value for the directional shear in deg/m
2.00                                    REFHEIGHT          - the reference height, used to contruct the BL profile

----------------------------------------Ocean Depth, Waves and Currents-------------------------------------------
the following parameters only need to be set if ISOFFSHORE = 1
200.00                                   WATERDEPTH         - the water depth
New_Wave.lwa                             WAVEFILE           - the path to the wave file, leave blank if unused
1                                        WAVESTRETCHING     - the type of wavestretching, 0 = vertical, 1 = wheeler, 2 = extrapolation, 3 = none
10000.00                                 SEABEDSTIFF        - the vertical seabed stiffness [N/m^3]
0.20                                     SEABEDDAMP         - a damping factor for the vertical seabed stiffness evaluation, between 0 and 1 [-]
0.10                                     SEABEDSHEAR        - a factor for the evaluation of shear forces (friction), between 0 and 1 [-]
0.00                                     SURF_CURR_U        - near surface current velocity [m/s]
0.00                                     SURF_CURR_DIR      - near surface current direction [deg]
30.00                                    SURF_CURR_DEPTH    - near surface current depth [m]
0.00                                     SUB_CURR_U         - sub surface current velocity [m/s]
0.00                                     SUB_CURR_DIR       - sub surface current direction [deg]
0.14                                     SUB_CURR_EXP       - sub surface current exponent
0.00                                     SHORE_CURR_U       - near shore (constant) current velocity [m/s]
0.00                                     SHORE_CURR_DIR     - near shore (constant) current direction [deg]

----------------------------------------Global Mooring System------------------------------------------------------
                                         MOORINGSYSTEM      - the path to the global mooring system file, leave blank if unused

----------------------------------------Dynamic Wake Meandering----------------------------------------------------
2                                        DWMSUMTYPE         - the dynamic wake meandering wake summation type: 0 = DOMINANT; 1 = QUADRATIC; 2 = LINEAR

----------------------------------------Environmental Parameters---------------------------------------------------
1.22500                                  DENSITYAIR         - the air density [kg/m^3]
0.000014607                              VISCOSITYAIR       - the air kinematic viscosity
1025.00000                               DENSITYWATER       - the water density [kg/m^3]
0.000001007                              VISCOSITYWATER     - the water kinematic viscosity [m^2/s]
9.806600000                              GRAVITY            - the gravity constant [m/s0000
----------------------------------------Output Parameters------------------------------0000------------------------
20.00000                                 STOREFROM          - the simulation stores data from this point in time, in [s]
false                                    STOREREPLAY        - store a replay of the simulation (warning, large memory will be required) [bool]
true                                     STOREAERO          - should the aerodynamic data be stored [bool]
false                                    STOREBLADE         - should the local aerodynamic blade data be stored [bool]
true                                     STORESTRUCT        - should the structural data be stored [bool]
true                                     STORESIM           - should the simulation (performance) data be stored [bool]
true                                     STOREHYDRO         - should the controller data be stored [bool]
false                                    STORECONTROLLER    - should the controller data be stored [bool]
false                                    STOREDWM           - should the dynamic wake meandering (DWM) data be stored [bool]

----------------------------------------Modal Analysis Parameters--------------------------------------------------
false                                    CALCMODAL          - perform a modal analysis (only single turbine simulations) [bool]
0.00000                                  MINFREQ            - store Eigenvalues, starting with this frequency
0.00000                                  DELTAFREQ          - omit Eigenvalues that are closer spaced than this value
100000000                                NUMFREQ            - set the number of Eigenmodes and Eigenvalues that will be stored"""
    simGen = open("./Generated/sim.txt", "w")
    simGen.writelines(sim)
    simGen.close()

def createolszakSim(windSpeed, rpm):
    sim = f"""----------------------------------------QBlade Simulation Definition File------------------------------------------
Generated with : QBlade CE v2.0.7-release_candidate_beta windows
Archive Format: 310023
Time : 21:27:48
Date : 15.05.2024

----------------------------------------Object Name-----------------------------------------------------------------
generated                            OBJECTNAME         - the name of the simulation object

----------------------------------------Simulation Type-------------------------------------------------------------
0                                        ISOFFSHORE         - use a number: 0 = onshore; 1 = offshore

----------------------------------------Turbine Parameters---------------------------------------------------------
multiple turbines can be added by adding multiple definitions encapsulated with TURB_X and END_TURB_X, where X must start at 1

TURB_1
    olszaktrb.txt                     TURBFILE           - the turbine definition file that is used for this simulation
    NREL_5MW_OC4_SEMI_RWT                TURBNAME           - the (unique) name of the turbine in the simulation (results will appear under this name)
    0                                INITIAL_YAW        - the initial turbine yaw in [deg]
    0                                INITIAL_PITCH      - the initial collective blade pitch in [deg]
    0                                INITIAL_AZIMUTH    - the initial azimuthal rotor angle in [deg]
    1                                    STRSUBSTEP         - the number of structural substeps per timestep (usually 1)
    1                                    RELAXSTEPS         - the number of initial static structural relaxation steps
    1                                    PRESCRIBETYPE      - rotor RPM prescribe type (0 = ramp-up; 1 = whole sim; 2 = no RPM prescibed)
    {rpm}                                RPMPRESCRIBED      - the prescribed rotor RPM [-]
    5                                    STRITERATIONS      - number of iterations for the time integration (used when integrator is HHT or Euler)
    0                                    MODNEWTONITER      - use the modified newton iteration?
    1                                    INCLUDEAERO        - include aerodynamic forces?
    1                                    INCLUDEHYDRO       - include hydrodynamic forces?
    0.00                                 GLOBPOS_X          - the global x-position of the turbine [m]
    0.00                                 GLOBPOS_Y          - the global y-position of the turbine [m]
    0.00                                 GLOBPOS_Z          - the global z-position of the turbine [m]
    0.00                                 GLOBROT_X          - the global x-rotation of the turbine [deg]
    0.00                                 GLOBROT_Y          - the global y-rotation of the turbine [deg]
    0.00                                 GLOBROT_Z          - the global z-rotation of the turbine [deg]
                                         EVENTFILE          - the file containing fault event definitions (leave blank if unused)
                                         LOADINGFILE        - the loading file name (leave blank if unused)
                                         SIMFILE            - the simulation file name (leave blank if unused)
                                         MOTIONFILE         - the prescribed motion file name (leave blank if unused)
END_TURB_1

----------------------------------------Simulation Settings-------------------------------------------------------
0.00050000                                 TIMESTEP           - the timestep size in [s]
200                                        NUMTIMESTEPS       - the number of timesteps
0.000                                    RAMPUP             - the rampup time for the structural model
0.000                                    ADDDAMP            - the initial time with additional damping
50.000                                  ADDDAMPFACTOR      - for the additional damping time this factor is used to increase the damping of all components
0.000                                    WAKEINTERACTION    - in case of multi-turbine simulation the wake interaction start at? [s]

----------------------------------------Wind Input-----------------------------------------------------------------
0                                        WNDTYPE            - use a number: 0 = steady; 1 = windfield; 2 = hubheight
                                         WNDNAME            - filename of the turbsim input file, mann input file or hubheight file (with extension), leave blank if unused
0                                        STITCHINGTYPE      - the windfield stitching type; 0 = periodic; 1 = mirror
true                                     WINDAUTOSHIFT      - the windfield shifting automatically based on rotor diameter [bool]
0.00                                     SHIFTTIME          - the windfield is shifted by this time if WINDAUTOSHIFT = 0
{windSpeed}                                    MEANINF            - the mean inflow velocity, overridden if a windfield or hubheight file is use
0.00                                     HORANGLE           - the horizontal inflow angle
0.00                                     VERTANGLE          - the vertical inflow angle
0                                        PROFILETYPE        - the type of wind profile used (0 = Power Law; 1 = Logarithmic)
0.00                                    SHEAREXP           - the shear exponent if using a power law profile, if a windfield is used these values are used to calculate the mean wake convection velocities
0.010                                    ROUGHLENGTH        - the roughness length if using a log profile, if a windfield is used these values are used to calculate the mean wake convection velocities
0.00                                     DIRSHEAR           - a value for the directional shear in deg/m
2.00                                    REFHEIGHT          - the reference height, used to contruct the BL profile

----------------------------------------Ocean Depth, Waves and Currents-------------------------------------------
the following parameters only need to be set if ISOFFSHORE = 1
200.00                                   WATERDEPTH         - the water depth
New_Wave.lwa                             WAVEFILE           - the path to the wave file, leave blank if unused
1                                        WAVESTRETCHING     - the type of wavestretching, 0 = vertical, 1 = wheeler, 2 = extrapolation, 3 = none
10000.00                                 SEABEDSTIFF        - the vertical seabed stiffness [N/m^3]
0.20                                     SEABEDDAMP         - a damping factor for the vertical seabed stiffness evaluation, between 0 and 1 [-]
0.10                                     SEABEDSHEAR        - a factor for the evaluation of shear forces (friction), between 0 and 1 [-]
0.00                                     SURF_CURR_U        - near surface current velocity [m/s]
0.00                                     SURF_CURR_DIR      - near surface current direction [deg]
30.00                                    SURF_CURR_DEPTH    - near surface current depth [m]
0.00                                     SUB_CURR_U         - sub surface current velocity [m/s]
0.00                                     SUB_CURR_DIR       - sub surface current direction [deg]
0.14                                     SUB_CURR_EXP       - sub surface current exponent
0.00                                     SHORE_CURR_U       - near shore (constant) current velocity [m/s]
0.00                                     SHORE_CURR_DIR     - near shore (constant) current direction [deg]

----------------------------------------Global Mooring System------------------------------------------------------
                                         MOORINGSYSTEM      - the path to the global mooring system file, leave blank if unused

----------------------------------------Dynamic Wake Meandering----------------------------------------------------
2                                        DWMSUMTYPE         - the dynamic wake meandering wake summation type: 0 = DOMINANT; 1 = QUADRATIC; 2 = LINEAR

----------------------------------------Environmental Parameters---------------------------------------------------
1.22500                                  DENSITYAIR         - the air density [kg/m^3]
0.000014607                              VISCOSITYAIR       - the air kinematic viscosity
1025.00000                               DENSITYWATER       - the water density [kg/m^3]
0.000001007                              VISCOSITYWATER     - the water kinematic viscosity [m^2/s]
9.806600000                              GRAVITY            - the gravity constant [m/s0000
----------------------------------------Output Parameters------------------------------0000------------------------
20.00000                                 STOREFROM          - the simulation stores data from this point in time, in [s]
false                                    STOREREPLAY        - store a replay of the simulation (warning, large memory will be required) [bool]
true                                     STOREAERO          - should the aerodynamic data be stored [bool]
false                                    STOREBLADE         - should the local aerodynamic blade data be stored [bool]
true                                     STORESTRUCT        - should the structural data be stored [bool]
true                                     STORESIM           - should the simulation (performance) data be stored [bool]
true                                     STOREHYDRO         - should the controller data be stored [bool]
false                                    STORECONTROLLER    - should the controller data be stored [bool]
false                                    STOREDWM           - should the dynamic wake meandering (DWM) data be stored [bool]

----------------------------------------Modal Analysis Parameters--------------------------------------------------
false                                    CALCMODAL          - perform a modal analysis (only single turbine simulations) [bool]
0.00000                                  MINFREQ            - store Eigenvalues, starting with this frequency
0.00000                                  DELTAFREQ          - omit Eigenvalues that are closer spaced than this value
100000000                                NUMFREQ            - set the number of Eigenmodes and Eigenvalues that will be stored"""
    simGen = open("./Generated/olszaksim.txt", "w")
    simGen.writelines(sim)
    simGen.close()

trbGen = open("./Generated/trb.txt", "w")
trbGen.writelines(trb)
trbGen.close()

createBlade([0.1, 0.074, 0.053, 0.042, 0.034, 0.029, 0.025, 0.022, 0.02, 0.018], [32.5, 18.6, 12.1, 8.4, 6, 4.4, 3.2, 2.2, 1.5, 0.9], "S1223", 120, 2)
#createSim(10, 0.5)
