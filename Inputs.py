import re

import numpy as np
from array import array
import parse
from parse import compile


class Inputs:
    def __init__(self, file_address):
        self.file = open(file_address, 'r').read()
        self.demands = []
        self.H0 = 0
        self.lengths = []
        self.diameters = []
        self.roughness = 0
        self.pipe_num = 0
        self.node_num = 0
        self.A = []


    def find_and_map(self):
        # Demands
        demands = self.file.split('\n\n')[1].split('\n')[3: len(self.file.split('\n\n')[1].split('\n'))]
        for i in range(len(demands)):
            self.demands.append(float(re.split(r'\s', demands[i].replace(' ', ''))[2]))
        self.node_num = len(demands) + 1
        print(self.node_num)

        # H0
        self.H0 = self.file.split('\n\n')[2].split('\n')[2].replace(' ', '').split('\t')[1]

        # Diameters - Lengths - Roughness
        x = self.file.split('\n\n')[4].split('\n')[2: len(self.file.split('\n\n')[4].split('\n'))]
        for i in range(len(x)):
            self.lengths.append(float(re.split(r'\s', x[i].replace(' ', ''))[3]))
            self.diameters.append(float(re.split(r'\s', x[i].replace(' ', ''))[4]))
        self.roughness = re.split(r'\s', x[0].replace(' ', ''))[5]
        print(self.lengths)
        self.pipe_num = len(x)
        print(self.pipe_num)


    def make_A12(self):
        edges = []
        x = self.file.split('\n\n')[4].split('\n')[2: len(self.file.split('\n\n')[4].split('\n'))]
        for i in range(len(x)):
            edges.append(
                (int(re.split(r'\s', x[i].replace(' ', ''))[1]), int(re.split(r'\s', x[i].replace(' ', ''))[2])))
        print(edges)
        for i in range(len(edges)):
            x = np.zeros([self.node_num])
            x[edges[i][0] - 1] = -1
            x[edges[i][1] - 1] = 1
            self.A.append(x)
        self.A = np.array(self.A)
        print(self.A)
