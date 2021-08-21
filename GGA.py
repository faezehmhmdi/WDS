import numpy as np
from Partition import Partition
import re

np.set_printoptions(suppress=True)

"""
The GGA Class
"""


class GGA:
    def __init__(self, file_address):
        # file
        self.file = open(file_address, 'r').read()

        # node_index Dictionary
        self.node_index = dict()

        # H0
        self.H0 = float(self.file.split('\n\n')[2].split('\n')[2].replace(' ', '').split('\t')[1])
        self.node_index[self.file.split('\n\n')[2].split('\n')[2].replace(' ', '').split('\t')[0]] = len(
            self.node_index)
        self.H0 = np.array(self.H0).reshape(-1, 1)

        # Demands (q)
        self.q_star = []
        demands = self.file.split('\n\n')[1].split('\n')[3: len(self.file.split('\n\n')[1].split('\n'))]
        for i in range(len(demands)):
            self.q_star.append(float(re.split(r'\s', demands[i].replace(' ', ''))[2]) / 1000)
            self.node_index[re.split(r'\s', demands[i].replace(' ', ''))[0]] = len(self.node_index)

        # D - C - L
        self.pipe_length = []
        self.pipe_diameter = []
        x = self.file.split('\n\n')[4].split('\n')[2: len(self.file.split('\n\n')[4].split('\n'))]
        for i in range(len(x)):
            self.pipe_length.append(float(re.split(r'\s', x[i].replace(' ', ''))[3]))
            self.pipe_diameter.append(float(re.split(r'\s', x[i].replace(' ', ''))[4]) / 1000)
        self.pipe_roughness = float(re.split(r'\s', x[0].replace(' ', ''))[5])

        # n
        self.n = 1.852

        # Np
        self.pipe_num = len(x)

        # Nn
        self.node_num = len(demands) + 1

        # A Matrix
        self.A = []
        self.make_A()

        # A10 Matrix
        self.A10 = np.array(self.A[:, 0]).reshape(-1, 1)

        # A12 Matrix
        self.A12 = np.array(self.A[:, 1:])

        # A21 Matrix
        self.A21 = self.A12.T

        # A22 Matrix
        self.A22 = np.zeros([self.node_num - 1, self.node_num - 1])

        # K
        self.K = []
        for i in range(self.pipe_num):
            k = (10.67 * self.pipe_length[i]) / (
                    self.pipe_roughness ** self.n * self.pipe_diameter[i] ** 4.87)
            self.K.append(float(k))

        # Q
        self.Q = self.calculate_Q()

        # D11
        self.D11 = np.diag(self.K * np.abs(self.Q) ** (self.n - 1))

        # H
        self.H = self.calculate_H()

        # F
        self.F = np.zeros([self.pipe_num + self.node_num - 1, 1])

        # b
        self.b = np.zeros([self.pipe_num + self.node_num - 1, 1])

        self.q_h = np.concatenate([np.array(np.abs(self.Q)).reshape(-1, 1), np.array(np.abs(self.H)).reshape(-1, 1)])

        # DA
        self.DA = np.zeros([self.node_num + self.pipe_num - 1, self.node_num + self.pipe_num - 1])

        # delta_QH
        self.deltaQ_deltaH = np.zeros([self.node_num + self.pipe_num - 1, 1])

        # Main Matrix
        self.main_matrix = np.block([[self.D11, self.A12], [self.A21, self.A22]])

        # Partition
        self.partition = Partition()

        # Subsystems
        self.subsystems = self.partition.subsystems

    def make_A(self):
        edges = []
        temp = self.file.split('\n\n')[4].split('\n')[2: len(self.file.split('\n\n')[4].split('\n'))]
        for i in range(len(temp)):
            edges.append(
                (self.node_index.get(re.split(r'\s', temp[i].replace(' ', ''))[1]),
                 self.node_index.get(re.split(r'\s', temp[i].replace(' ', ''))[2])))
        for i in range(len(edges)):
            temp = np.zeros([self.node_num])
            temp[edges[i][0]] = -1
            temp[edges[i][1]] = 1
            self.A.append(temp)
        self.A = np.array(self.A)

    def calculate_Q(self):
        x = []
        a = self.A21
        b = np.array(self.q_star).reshape(-1, 1)
        for i in range(len(np.linalg.lstsq(a, b, rcond=None)[0])):
            x.append(np.linalg.lstsq(a, b, rcond=None)[0][i][0])
        return x

    def calculate_H(self):
        x = []
        a = self.A12
        b = -self.A10 @ self.H0 - self.D11 @ np.array(self.Q).reshape(-1, 1)
        for i in range(len(np.linalg.lstsq(a, b, rcond=None)[0])):
            x.append(np.linalg.lstsq(a, b, rcond=None)[0][i][0])
        return x

    def gga_algorithm(self):
        count = 0
        self.b = np.concatenate([-self.A10 @ self.H0, np.array(self.q_star).reshape(-1, 1)])
        self.F = np.concatenate(
            [self.D11 @ np.array(np.abs(self.Q)).reshape(-1, 1) + self.A12 @ np.array(np.abs(self.H)).reshape(-1, 1),
             self.A21 @ np.array(np.abs(self.Q)).reshape(-1, 1)]) - self.b
        # print('The Initial F Matrix: \n')
        # print(self.F)
        dE = np.ones([self.node_num, 1])
        dq = np.ones([self.node_num, 1])
        while abs(np.linalg.norm(dE)) >= 0.01:
            self.D11 = self.n * self.D11
            self.DA = np.block([[self.D11, self.A12], [self.A21, self.A22]])
            self.deltaQ_deltaH = -(np.linalg.inv(self.DA)) @ self.F
            print('\nThe deltaQ deltaH Matrix in iteration: ' + str(count) + '\n')
            print(self.deltaQ_deltaH)
            self.q_h += self.deltaQ_deltaH
            print('\nThe updated QH Matrix in iteration: ' + str(count) + '\n')
            print(self.q_h)
            # print(self.q_h[0:self.pipe_num].T[0])
            self.D11 = np.diag(self.K * np.abs(self.q_h[0:self.pipe_num]).T[0] ** (self.n - 1))
            self.q_star = self.A21 @ self.q_h[0:self.pipe_num]
            self.b = np.concatenate([-self.A10 @ self.H0, self.q_star])
            self.F = np.concatenate(
                [self.D11 @ self.q_h[0:self.pipe_num] + self.A12 @ self.q_h[
                                                                   self.pipe_num:self.pipe_num + self.node_num],
                 self.A21 @ self.q_h[0:self.pipe_num]]) - self.b
            dE = -self.F[0:self.pipe_num]
            print('\nThe updated dE Matrix in iteration: ' + str(count) + '\n')
            print(dE)
            dq = -self.F[self.pipe_num:self.pipe_num + self.node_num]
            # print('\nThe updated dq Matrix in iteration: ' + str(count) + '\n')
            # print(dq)
            count += 1
        # print(count)
