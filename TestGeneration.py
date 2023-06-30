import numpy as np

row = 1
col = 1000

x = []
x = np.linspace(0, 6, col).tolist()
f = open("./Output/Test-Input.txt", "w")
f.write(str(row)+" ")
f.write(str(col)+" ")
for i in x:
    f.write(str(i)+" ")
f.close()