import matplotlib.pyplot as plt
from datetime import datetime

#cut=False
cut=True

plt.figure(figsize=(16, 10))
#plt.subplots_adjust(wspace = 0.2)
plt.subplots_adjust(hspace = 0.3)

t1, temp = [], []
for line in open('res.dat', 'r'):
    lines = [float(s) for s in line.split()]
    dt_object = datetime.fromtimestamp(lines[0])
    #print("dt_object =", dt_object)
     
    if cut==True and ( ( lines[0] > 1716366300 and lines[0] < 1716370500 ) or ( lines[0] > 1716372900 and lines[0] < 1716373200 )
            or ( lines[0] > 1716433200 and lines[0] < 1716442800 ) or ( lines[0] > 1716453000 and lines[0] < 1716460200 )
            or ( lines[0] > 1716514800 and lines[0] < 1716525720 ) or ( lines[0] > 1716526500 and lines[0] < 1716532200 )
            or ( lines[0] > 1716539400 and lines[0] < 1716545280 ) or ( lines[0] > 1716545460 and lines[0] < 1716546600 ) 
            or ( lines[0] > 1716780600 and lines[0] < 1716787800 ) or ( lines[0] > 1716793680 and lines[0] < 1716801900 ) 
            or ( lines[0] > 1716804300 and lines[0] < 1716805200 ) or ( lines[0] > 1716865980 and lines[0] < 1716876300 )
            or ( lines[0] > 1716878100 and lines[0] < 1716897600 ) or ( lines[0] > 1716947400 and lines[0] < 1716949800 )
            or ( lines[0] > 1716954000 and lines[0] < 1716962100 ) or ( lines[0] > 1716964500 and lines[0] < 1716978300 ) 
            or ( lines[0] > 1717037459 and lines[0] < 1717040459 ) or ( lines[0] > 1717046100 and lines[0] < 1717049580 ) 
            or ( lines[0] > 1717050900 and lines[0] < 1717393440) or ( lines[0] > 1717394880 ) ):
        t1.append(dt_object)
        temp.append(lines[1])
    elif cut==False:
        t1.append(dt_object)
        temp.append(lines[1])

listFiles = ["channel_0.dat","channel_2.dat","channel_4.dat"]
#listFiles = ["channel_0.dat","channel_2.dat","channel_4.dat","noise_channel_0.dat"]

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
        
        if cut==True and ( ( lines[0] > 1716366300 and lines[0] < 1716370500 ) or ( lines[0] > 1716372900 and lines[0] < 1716373200 ) 
                or ( lines[0] > 1716433200 and lines[0] < 1716442800 ) or ( lines[0] > 1716453000 and lines[0] < 1716460200 )
                or ( lines[0] > 1716514800 and lines[0] < 1716525720 ) or ( lines[0] > 1716526500 and lines[0] < 1716532200 )
                or ( lines[0] > 1716539400 and lines[0] < 1716545280 ) or ( lines[0] > 1716545460 and lines[0] < 1716546600 )
                or ( lines[0] > 1716780600 and lines[0] < 1716787800 ) or ( lines[0] > 1716793680 and lines[0] < 1716801900 )
                or ( lines[0] > 1716804300 and lines[0] < 1716805200 ) or ( lines[0] > 1716865980 and lines[0] < 1716876300 ) 
                or ( lines[0] > 1716878100 and lines[0] < 1716897600 ) or ( lines[0] > 1716947400 and lines[0] < 1716949800 )
                or ( lines[0] > 1716954000 and lines[0] < 1716962100 ) or ( lines[0] > 1716964500 and lines[0] < 1716978300 )
                or ( lines[0] > 1717037459 and lines[0] < 1717040459 ) or ( lines[0] > 1717046100 and lines[0] < 1717049580 )
                or ( lines[0] > 1717050900 and lines[0] < 1717393440) or ( lines[0] > 1717394880 ) ):
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
plt.axvline(x=datetime.fromtimestamp(1716804600), label='to light guide', c='g')
plt.axvline(x=datetime.fromtimestamp(1717043400), label='thermofuse', c='b')
plt.axvline(x=datetime.fromtimestamp(1717393800), label='set up isotope Co', c='r')
plt.xticks(rotation = 25)
#plt.autoscale(enable=True, axis='both', tight = None)

plt.subplot(212)
#plt.title("....")
plt.xlabel('Time')
plt.ylabel('LED amplitude')
plt.legend()
plt.ylim(0, 100)
plt.axvline(x=datetime.fromtimestamp(1716804600), label='to light guide', c='g')
plt.axvline(x=datetime.fromtimestamp(1717043400), label='thermofuse', c='b')
plt.axvline(x=datetime.fromtimestamp(1717393800), label='set up isotope Co', c='r')
plt.xticks(rotation = 25)
#plt.autoscale(enable=True, axis='both', tight = None)

plt.show()

