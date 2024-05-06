import numpy as np
from functools import reduce
import time


def readxyz(filename, times, molecules):
    file1 = open(filename)
    finxyz = [[[] for _ in range(molecules)] for _ in range(times)]
    for pertime in range(times):
        t1 = file1.readline()
        for moles in range(molecules):
            perxyz = [float(ix) for ix in file1.readline().split()[-3:]]
            finxyz[pertime][moles] = [perxyz[0],perxyz[2],perxyz[1]]      #从前到后为xyz,算的2对应的方向
    file1.close()
    return np.array(finxyz, dtype=np.float32)


def xshift(watt1, watt2):
    return (watt1[0] - watt2[0])*(watt1[0] - watt2[0])


def yshift(watt1, watt2):
    return (watt1[1] - watt2[1])*(watt1[1] - watt2[1])


def msdxy(max_time, re_left, re_right, arr_multi):
    fin_dip_x = []
    fin_dip_y = []
    fin_dip_down = []
    datashape = arr_multi.shape
    per_pair = [[] for _ in range(datashape[1])]
    bin_in = np.zeros((datashape[0], datashape[1]))
    for i in range(datashape[0]):
        for j in range(datashape[1]):
            if re_left < arr_multi[i][j][2] < re_right:
                bin_in[i][j] = 1
                per_pair[j].append([i, i])
    fin_dip_x.append(0)
    fin_dip_y.append(0)
    fin_dip_down.append(np.sum(bin_in))
    for time in range(1, min([int(datashape[0] / 2), max_time])):
        diplist_x = []
        diplist_y = []
        diplist_2 = []
        for wat_ind in range(datashape[1]):
            remove_list = []
            for per_pre in range(len(per_pair[wat_ind])):
                ftime = per_pair[wat_ind][per_pre][0]
                nowtime = per_pair[wat_ind][per_pre][1] + 1
                if nowtime < datashape[0] and bin_in[nowtime][wat_ind]:
                    diplist_x.append(xshift(watt1=arr_multi[ftime][wat_ind], watt2=arr_multi[nowtime][wat_ind]))
                    diplist_y.append(yshift(watt1=arr_multi[ftime][wat_ind], watt2=arr_multi[nowtime][wat_ind]))
                    diplist_2.append(1)
                    per_pair[wat_ind][per_pre][1] += 1
                else:
                    remove_list.append(per_pair[wat_ind][per_pre])
            for i in remove_list:
                per_pair[wat_ind].remove(i)
        if len(diplist_x):
            fin_dip_x.append(reduce(lambda x, y: x + y, diplist_x))
            fin_dip_y.append(reduce(lambda x, y: x + y, diplist_y))
            fin_dip_down.append(reduce(lambda x, y: x + y, diplist_2))
    return fin_dip_x, fin_dip_y, fin_dip_down


def cut_piece(kk, piece):
    fin = []
    while kk[1] - kk[0] > piece*2:
        fin.append([kk[0], kk[0]+piece])
        fin.append([kk[1]-piece, kk[1]])
        kk[0] += piece
        kk[1] -= piece
    fin.append([kk[0], kk[1]])
    fin.sort()
    return fin


def msdout(max_time, re_left, re_right, arr_multi):
    namelat = "_" + ('000' + str(int(re_left)))[-4:] + "_" + ('000' + str(int(re_right)))[-4:]
    fname = open("msd" + namelat, "w")
    fgh = msdxy(max_time, re_left, re_right, arr_multi)
    for io in range(len(fgh[0])):
        print("%.5f %.5f %.5f " % (fgh[0][io], fgh[1][io], fgh[2][io]), sep=" ", file=fname)
    fname.close()


def finmsd():
    piece1 = 90 #bin间隔,单位A
    #cl_left, cl_right = 49.13, 170.83
    cl_left, cl_right = 0, 90 #上下总边界
    totalmols = 8212  #总分子数
    times = 4000  #帧数
    allmols = np.array([1, 1007, 635, 1000, 1000, 1000, 1000, 1000, 1569])  #各类分子数
    # allmols = np.array([1, 5000, 95, 188, 142, 113, 142, 113, 220, 306, 273, 151, 48, 32, 90, 47, 47, 43])
    mol_index = 1 # 1是算array中第二个，2是算第三个 得出第一二行除第三为x y msd

    mol_left = np.sum(allmols[:mol_index]) + 1
    mol_right = np.sum(allmols[:mol_index+1])
    print(mol_left, mol_right)
    kk1 = [cl_left, cl_right]
    fin1 = cut_piece(kk1, piece1)
    arr_multi = readxyz('newcen', 4000, 8212) #名字,帧数,分子数
    for ix in range(len(fin1)):
        msdout(1500, fin1[ix][0], fin1[ix][1], arr_multi[:, mol_left-1:mol_right, :]) #out-msdtime


t1 = time.time()
finmsd()
t2 = time.time()
print(t2-t1)
