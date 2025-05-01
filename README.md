# Individual-Project

This project is used to calculate the optimal blade design for use in a mechanical transmission wind powered vehicle with the purpose of achieveing maximum Wach number. It is developed and used within the individual project conducted by Alex Vale at the University of Southampton. This contains a full documentation of the code and how it works as well as an application from a starting blade.

Running the Program:
1. Make sure QBlade is installed on the device in the executable folder. Visit https://qblade.org/ for more help.
2. For your chosen airfoil, use QBlade to get 360 polars for a range of reynolds numbers (typically 50k to 1000k) and store these in a folder local to the executable.
3. Ensure physical constants in "race-funct.py" match those of your chosen vehicle.
4. Ensure the name of your selected airfoil is correctly inputed in the "fullcode.py" script.
5. Input starting parameters in the "x0" variable in the "fullcode.py" script. The first 20 describe blade geometry inputs (10 chord then 10 twist) at 10 equally spaced positions along the blade. The final 5 inputs are initial conditions for the gearbox change. The first 2 are the velocities at which gear changes occur and the final 3 are the gears in order.
6. Make a folder local to the executable named "logs".
7. The code can now be run using the "fullcode.py" script. WARNING: THIS MAY TAKE SOME TIME DEPENDING ON DEVICE CAPABILITY

Reading Data:
Each executed optimisation is allocated a random 6 digit hex code to label it. For data extraction, make a note of this from the logs folder and input into the solve variable in the following scripts.

The "Extract Optimal Blade Description.py" script can be run to generate a QBlade .bld file which describes the optimal blade design. This can be imported into QBlade to visualise the blade.

"Log Reader 2.0" allows you to look at the data for any particular blade iteration within the solve. Adapt the script as necessary to get the output graphs desired.

"Control Points log" Tracks the movement of the objective variables as the optimiser improves the design.
