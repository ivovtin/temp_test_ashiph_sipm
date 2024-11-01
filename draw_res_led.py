import matplotlib.pyplot as plt
from datetime import datetime

plt.figure(figsize=(16, 10))
plt.subplots_adjust(wspace = 0.2)

cut=False
#cut=True

t1, temp = [], []
for line in open('res.dat', 'r'):
    lines = [float(s) for s in line.split()]
    dt_object = datetime.fromtimestamp(lines[0])
    #print("dt_object =", dt_object)
     
    if cut==True and lines[0] > 1716779197:
        t1.append(dt_object)
        temp.append(lines[1])
    elif cut==False:
        t1.append(dt_object)
        temp.append(lines[1])

listFiles = ["channel_0.dat","channel_2.dat","channel_4.dat"]

for fileName in listFiles:
    t2, Amp, errAmp = [], [], []
    for line in open(fileName, 'r'):
        lines = [float(s) for s in line.split()]
        dt_object = datetime.fromtimestamp(lines[0])

        A1phe = []
        if listFiles.index(fileName) == 0: 
            A1phe = 0.8*pow(10,-6)
            ilabel = 'Channel 0'
        elif listFiles.index(fileName) == 1:
            A1phe = 0.5*pow(10,-6)
            if lines[0] > 1716512830:
                A1phe = 0.8*pow(10,-6)
            ilabel = 'Channel 2'
        elif listFiles.index(fileName) == 2:
            A1phe = 0.35*pow(10,-6)
            if lines[0] > 1716512830:
                A1phe = 0.8*pow(10,-6)
            if lines[0] > 1716791640 and lines[0] < 1716803700:
                A1phe = 0.5*pow(10,-6)
            if lines[0] > 1717316741:
                A1phe = 0.95*pow(10,-6)
            ilabel = 'Channel 4'
        elif listFiles.index(fileName) == 3:
            A1phe = 0.8*pow(10,-6)
            ilabel = 'Ped Ch0'
        
        if cut==True and lines[0] > 1716779197:
            t2.append(dt_object)
            Amp.append(lines[1]/A1phe)
            errAmp.append(lines[2]/(10000**0.5*A1phe))
        elif cut==False:
            t2.append(dt_object)
            Amp.append(lines[1]/A1phe)
            errAmp.append(lines[2]/(10000**0.5*A1phe))

    plt.subplot(212)
    plt.errorbar(t2, Amp, yerr=errAmp, marker = 'o', markersize=2, linestyle='None', label=ilabel)

plt.subplot(211)
#plt.title("....")
plt.xlabel('Time')
plt.ylabel('Temperature, °C')
plt.plot(t1, temp, marker = 'o', markersize=2, c = 'b', linestyle='None', linewidth=2)
plt.xticks(rotation = 25)
#plt.autoscale(enable=True, axis='both', tight = None)

plt.subplot(212)
#plt.title("....")
plt.xlabel('Time')
plt.ylabel('LED amplitude')
plt.legend()
plt.ylim(0, 100)
plt.xticks(rotation = 25)
#plt.autoscale(enable=True, axis='both', tight = None)

plt.show()

