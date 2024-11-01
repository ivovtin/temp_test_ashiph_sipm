import matplotlib.pyplot as plt
from datetime import datetime

t, temp = [], []
for line in open('res.dat', 'r'):
    lines = [float(s) for s in line.split()]
    dt_object = datetime.fromtimestamp(lines[0])
    #print("dt_object =", dt_object)
    t.append(dt_object)
    temp.append(lines[1])

#plt.title("....")
plt.xlabel('Time')
plt.ylabel('Temperature, °C')
plt.plot(t, temp, marker = 'o', markersize=2, c = 'b', linestyle='None', linewidth=2)
plt.xticks(rotation = 25)
plt.autoscale(enable=True, axis='both', tight = None)
plt.show()
