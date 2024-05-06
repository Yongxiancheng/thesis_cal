import numpy as np
from functools import reduce
import matplotlib.pyplot as plt
import time


def readcenxyz(fname):
    file1 = open(fname)
    frames = []
    while 1:
        newline = file1.readline()
        if newline:
            pframe = []
            atomnum = int(newline.split()[1])
            for ix in range(atomnum):
                pline = file1.readline().split()[1:]
                pframe.append([float(jx) for jx in pline])
            frames.append(pframe)
        else:
            return np.array(frames)


def xshift(watt1, watt2):
    return (watt1 - watt2)*(watt1 - watt2)


def msdxy(data1, max_time=1500): #迭代最大值
    fin_dip_up = []
    fin_dip_down = []
    for time in range(0, min([int(data1.shape[0] / 2) + 1, max_time])):
        fin_dip_1 = []
        fin_dip_2 = []
        for numb in range(data1.shape[0] - time):
            diplist_1 = []
            diplist_2 = []
            for wat_ind in range(data1.shape[1]):
                diplist_1.append(xshift(watt1=data1[numb][wat_ind], watt2=data1[numb + time][wat_ind]))
                diplist_2.append(1)
            fin_dip_1.append(reduce(lambda x, y: x + y, diplist_1))
            fin_dip_2.append(reduce(lambda x, y: x + y, diplist_2))
        fin_dip_up.append(reduce(lambda x, y: x + y, fin_dip_1))
        fin_dip_down.append(reduce(lambda x, y: x + y, fin_dip_2))
    return fin_dip_up, fin_dip_down


def getmsd(data1, tstep):
    xup, xdn = msdxy(data1, 1500) #迭代最大值
    finy = np.array([xup[ix]/xdn[ix] for ix in range(len(xup))])
    finx = np.array([ix*tstep for ix in range(len(xup))])
    resmsd = np.polyfit(finx, finy, 1)
    return resmsd[0], finy


t1 = time.time()
datas = readcenxyz("newcen")
timestep = 1 #1=1ps 每帧间隔时间长度
molist = [1, 985, 11, 11, 635, 985, 11, 11, 880, 10, 10, 4560, 51, 51]   #各类分子数
molsum = [sum(molist[:ix+1]) for ix in range(len(molist))]
molsum.append(0)
allmsd = []
sqr_dis = []
mass = [0, 18, 23, 35.5, 16.043, 18, 23, 35.5, 18, 23, 35.5, 18, 23, 35.5]  #对应分子数的质量
for ix in range(len(molsum)-1):  #for ix in range(5, len(molsum)-1):  从第五个分子后开始算
    msddata = getmsd(datas[:, molsum[ix-1]:molsum[ix], 2], timestep)  #计算各方向msd; 0=x,1=y,2=z
    allmsd.append(msddata[0])
    sqr_dis.append(msddata[1])
t2 = time.time()
print("time:", t2-t1)
fileout1 = open("allmsd", "w")
fileout2 = open("sqdis", "w")
for ix in range(len(mass)):                         #for ix in range(len(allmsd)):
    print(mass[ix], allmsd[ix], file=fileout1)      #    print(mass[5+ix], allmsd[ix], file=fileout1) 
pltx = [ix*timestep for ix in range(len(sqr_dis[0]))]
for ix in range(len(sqr_dis[0])):
    print(pltx[ix], end=' ', file=fileout2)
    for jx in range(len(sqr_dis)):
        print("%f" % sqr_dis[jx][ix], end=' ', file=fileout2)
    print("\n", end='', file=fileout2)
fileout1.close()
fileout2.close()

# plt.figure(1)
# plt.scatter(mass, allmsd)
# plt.figure(2)
# for ix in range(len(sqr_dis)):
#     plt.plot(pltx, sqr_dis[ix])
# plt.show()
# 算出的单方向msd系数未除2
