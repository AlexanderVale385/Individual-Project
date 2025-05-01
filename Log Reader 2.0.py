import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

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

solve = "a6aa76"

main = pd.read_csv(f'./logs/{solve}/mainBlade.csv', skiprows=0).to_numpy()
data = np.array([[manage(word) for word in row] for row in main]) #0: x, 1: wach, 2:error, 3: iteration

#looking for this condition
maxWach = np.max(data[:, 1])
for i in range(len(data[:, 0])):
    if data[i, 1] == maxWach:
        print(f"Max Wach = {maxWach}")
        iteration = data[i, 3]
        print(iteration)
        print(data[i, 0])

arrays = np.loadtxt(f"./logs/{solve}/1501.csv", delimiter= ",") #0: time, 1: vel, 2: tsr, 3:cq, 4:ct, 5:drag, 6:prop, 7:net, 8:dist
arrays2 = np.loadtxt(f"./logs/12b6ac/1099.csv", delimiter= ",")
arrays3 = np.loadtxt(f"./logs/olszak final 2/490.csv", delimiter= ",")
# plt.plot(arrays[8], arrays[7], label = "net force")
plt.plot(arrays3[8], arrays3[1], label = "Initial Blade")
plt.plot(arrays[8], arrays[1], label = "Original Method")
plt.plot(arrays2[8], arrays2[1], label = "Alternative Method")
# plt.plot(arrays[0], arrays[2], label = "tsr")
# plt.plot(arrays[0], -1 * arrays[5], label = "drag")
# plt.plot(arrays[0], arrays[6], label = "pro")
plt.axvline(x=100, color='r', linestyle='--', label='Race Start')
plt.ylabel("Velocity (m/s)")
plt.xlabel("Total Distance Travelled (m)")
plt.legend()
plt.savefig(f"./figures/alternative vs original velocity profile.png")
plt.show()


