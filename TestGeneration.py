import numpy as np

row = 1
col = 1000

x = []
x = np.linspace(-5, 5, col).tolist()
f = open("./Output/Test-Input.txt", "w")
f.write(str(row)+" ")
f.write(str(col)+" ")
for i in x[0:-1]:
    f.write(str(i)+" ")
f.write(str(x[-1]))
f.close()