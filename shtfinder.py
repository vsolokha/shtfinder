#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
module shtfinder for search in Globus-M plasma pulses
made for calculation the euclidian and DTW distance of the 
plasma discharge diagnostics time series and consecutive sorting

Created by V. Solokha
Last updated    2019-09-07: Function created
                2019-09-20: WINE used for linux usage of test.exe 
                2019-09-27: Improved input reading routine, merged with convert.py

==========
MODULES:
==========

pandas, numpy, matplotlib, scipy, fastdtw, tqdm

==========
USAGE:
==========
python3 shtfinder.py baseshot start_sequence end_sequence diagnostic_idx weights shtpath

Example: python3 36612 36610 36614 6,7 1.0,0.0 "./sht"

==========
INPUTS:
==========

baseshot (int)         -- the number of the base pulse
start_sequence (int)   -- starting point of the range
end_sequence (int)     -- finishing point of the range 
diagnostic_idx (array) -- indices of the comparable diagnostics 
weights (array)        -- weights of the diagnostics
shtpath (str)          -- path to the sht files folder

'''
import pandas as pd
import numpy as np
import matplotlib.pylab as plt
from scipy.spatial.distance import euclidean
from fastdtw import fastdtw
from io import StringIO
from tqdm import tqdm
import sys
import os

class Shot:
    def __init__(self, number, shtpath):
        self.number  = number
        self.shtpath = shtpath 
        self.read()

    def read(self):
        filename_sht = f"/sht{self.number}.SHT"
        filename = f"./sht{self.number}.csv"
        if os.name == 'nt':
            os.system(r'test.exe "{0}" "{1}"'.format(self.shtpath + filename_sht, filename))
        else:
            os.system(r'WINEDEBUG=fixme-all LANG=ru_RU.UTF-8 wine test.exe "{0}" "{1}"'.format(self.shtpath + filename_sht, filename))
        with open(filename, 'r', encoding = "ISO-8859-1") as fp:
            line = fp.readline()
            n_oscillo = int(line)
            self.names = []
            self.info  = []
            self.unit  = []
            self.data  = []
            for idx_oscillo in range(n_oscillo):
                line    = fp.readline()
                self.names.append(line)
                line    = fp.readline()
                self.info.append(line)
                line    = fp.readline()
                self.unit.append(line)
                line    = fp.readline()
                length  = int(line) 
                oscillo = np.zeros([2, length])  
                for idx in range(length):
                    line  = fp.readline().split()
                    for i in range(oscillo.shape[0]):
                        oscillo[i, idx] = float(line[i])            
                self.data.append(oscillo)
        os.remove(filename)
                
    def print_names(self):
        i = 0
        for name in self.names:
            print(i, name)
            i += 1

    def get_ip(self):
        ip_idx = 1
        return(np.sort(self.data[ip_idx][1])[-100:-1].mean())

    def plot(self, idx, fignum=None, color = 'k'):
        plt.figure(fignum)
        plt.grid(True)
        plt.xlabel("t, ms")
        plt.ylabel(self.unit[idx])
        plt.title(self.names[idx])
        plt.plot(self.data[idx][0], self.data[idx][1], color=color)

    def get_data(self, columns): 
        data = []
        for col in columns:
            try:
                data.append(self.data[col])
            except:
                data.append(np.zeros([9, 4096]))
        return data


def euc_dist(a,b):
    a = a - a.mean()
    b = b - b.mean()
    return euclidean(a,b)

def compare_shots(baseshot, testshot, columns, method):
    bdata = baseshot.get_data(columns)
    tdata = testshot.get_data(columns)
    dist = []
    for idx in range(len(bdata)):
        if method == "dtw":
            distance, _ = fastdtw(bdata[idx][1], tdata[idx][1], dist=euclidean)
        elif method == "euc": 
            try:
                distance = euc_dist(bdata[idx][1][:len(tdata[idx][1])], tdata[idx][1])
            except:
                distance = -1.0
        dist.append(distance)
    return dist

def get_distance(baseshot_number, testshot_range, columns, shtpath, method="dtw"):
    baseshot = Shot(baseshot_number, shtpath)
    number   = [] 
    result   = []
    for testshot_number in tqdm(range(testshot_range[0], testshot_range[1]+1)):
        try: 
            testshot = Shot(testshot_number, shtpath)
            dist = compare_shots(testshot, baseshot, columns, method)
            number.append(testshot_number)
            result.append(dist)
            del(testshot)
        except:
            pass
    return np.array(result), np.array(number)
    
def get_wsorted(number, res, weight):
    dist = np.zeros(res.shape)
    for idx in range(res.shape[0]):
        for idx_w in range(res.shape[1]):
            w = weight[idx_w]/np.max(weight)
            dist[idx, idx_w] = w * res[idx,idx_w]
    dist_flat = np.sum(dist, axis=1)
    return np.c_[number[np.argsort(dist_flat)], res[np.argsort(dist_flat)]]
        


def main():
    baseshot   = int(sys.argv[1])
    test_range = [int(sys.argv[2]), int(sys.argv[3])] 
    columns    = list(map(int,   sys.argv[4].split(",")))
    weights    = list(map(float, sys.argv[5].split(",")))
    shtpath    = str(sys.argv[6])
    print(sys.argv)
    print(f"BASESHOT: {baseshot}")
    print(f"TESTSHOT RANGE: {test_range}")
    res_dtw, __     = get_distance(baseshot, test_range, columns, shtpath, "dtw")
    res_euc, number = get_distance(baseshot, test_range, columns, shtpath, "euc")
    res_dtw_sort = get_wsorted(number, res_dtw, weights)
    np.savetxt(f"result-dtw-{test_range[0]}.txt", res_dtw_sort, \
        fmt=['%6d','%5.4e','%5.4e'], header='DTW DISTANCE') 
    res_euc_sort = get_wsorted(number, res_euc, weights)
    np.savetxt(f"result-euc-{test_range[0]}.txt", res_euc_sort, \
        fmt=['%6d','%5.4e','%5.4e'], header='EUCLIDIAN DISTANCE')


if __name__ == "__main__":
    main()
