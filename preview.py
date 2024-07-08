import matplotlib.pyplot as plt
import numpy as np

In  = './Output/Test-Input.txt'
Out = './Output/Test-Output.txt'

fp = open(In, "r")
x = []
for line in fp.readlines():
    line=line.replace('\n','')
    line=line.split(' ')
    x.append(line)
fp.close()

del x[0][0]
del x[0][0]
x = x[0]

fp = open(Out, "r")
y1 = y2 = []
for line in fp.readlines():
    line=line.replace('\n','')
    line=line.split(' ')
    y1.append(line)            
fp.close()

y1 = y1[0]

x = list(map(float,x))
x = [i+1 for i in x]

y1 = list(map(float,y1))
y2 = np.sin(x)

plt.plot(x, y1, color='red')
plt.plot(x, y2, color='skyblue')

# 右上方对曲线颜色解释红色为测试输出，蓝色为标准输出
plt.legend(['Prediction', 'sin(x+1)'], loc='upper right')
plt.xlabel('x')
plt.ylabel('y')

plt.show()