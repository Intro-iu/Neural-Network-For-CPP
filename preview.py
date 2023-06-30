import matplotlib.pyplot as plt
import numpy as np

In  = './Output/Test-Input.txt'
Out = './Output/Test-Output.txt'

fp = open(In, "r")
x = []
y3 = []
for line in fp.readlines():
    line=line.replace('\n','')
    line=line.split(' ')
    x.append(line)
fp.close()

del x[0][0]
del x[0][0]
del x[0][-1]
x = x[0]

fp = open(Out, "r")
y1 = y2 = []
for line in fp.readlines():
    line=line.replace('\n','')
    line=line.split(' ')
    y1.append(line)            
fp.close()

del y1[0][-1]
del y1[1]
y1 = y1[0]

x = list(map(float,x))
y1 = list(map(float,y1))
y2 = np.sin(x)+6
for i in x:
    y3.append(i*i)


plt.scatter(x, y1, marker='x', color='springgreen')
plt.plot(x, y2, color='whitesmoke')
plt.scatter(x, y3, marker='x', color='royalblue')
plt.show()