import numpy as np
from Partition import Partition


np.set_printoptions(suppress=True)

"""
The GGA Class
"""


class GGA:
    def __init__(self):
        # n
        self.n = 1.852
        # Np
        self.pipe_num = 7
        # Nn
        self.node_num = 6
        # L
        self.pipe_length = 1000
        # C
        self.pipe_roughness = 100
        # d of pipes
        self.pipe_diameter = [0.6, 0.6, 0.5, 0.4, 0.3, 0.3, 0.3]
        # Q
        self.Q = np.array(list(map(float, '0.17 0.17 0.10 0.17 0.05 0.08 0.05'.split())))
        # H
        self.H = np.array(list(map(float, '35.35 35.00 34.00 28.00 22.00 23.00'.split())))
        # H0 Vector
        self.H0 = np.array([self.H[0]]).reshape(-1, 1)
        # K
        self.K = np.array(list(map(float, '25.397 25.397 61.728 183.033 743.202 743.202 743.202'.split())))
        # F
        self.F = np.zeros([self.pipe_num + self.node_num - 1, 1])
        # b
        self.b = np.zeros([self.pipe_num + self.node_num - 1, 1])
        # A10 Matrix
        self.A10 = np.array([
            [-1],
            [-1],
            [0],
            [0],
            [0],
            [0],
            [0]
        ])
        # D11
        self.D11 = np.diag(self.K * self.Q ** (self.n - 1))
        # A12 Matrix
        self.A12 = np.matrix([[1, 0, 0, 0, 0],
                              [0, 1, 0, 0, 0],
                              [-1, 1, 0, 0, 0],
                              [0, -1, 1, 0, 0],
                              [0, 0, -1, 1, 0],
                              [0, 0, -1, 0, 1],
                              [0, 0, 0, 1, -1]])
        # A21 Matrix
        self.A21 = self.A12.T
        # A22 Matrix
        self.A22 = np.zeros([self.node_num - 1, self.node_num - 1])
        # q_star
        self.Q = self.Q.reshape(-1, 1)
        self.q_star = self.A21 @ self.Q
        # q_h
        self.H = self.H[1:].reshape(-1, 1)
        self.q_h = np.concatenate([self.Q, self.H])
        # DA
        self.DA = np.zeros([self.node_num + self.pipe_num - 1, self.node_num + self.pipe_num - 1])
        # delta_QH
        self.deltaQ_deltaH = np.zeros([self.node_num + self.pipe_num - 1, 1])

        # Main Matrix
        self.main_matrix = np.block([[self.D11, self.A12], [self.A21, self.A22]])
        # print(self.main_matrix)

        # Partition
        self.partition = Partition()

        # Subsystems
        self.subsystems = self.partition.subsystems

    def gga_algorithm(self):
        count = 0
        self.b = np.concatenate([-self.A10 @ self.H0, self.q_star])
        self.F = np.concatenate(
            [self.D11 @ self.Q.reshape(-1, 1) + self.A12 @ self.H, self.A21 @ self.Q.reshape(-1, 1)]) - self.b
        # print('The Initial F Matrix: \n')
        # print(self.F)
        dE = np.ones([self.node_num, 1])
        dq = np.ones([self.node_num, 1])
        while abs(np.linalg.norm(dE)) >= 0.01:
            self.D11 = self.n * self.D11
            self.DA = np.block([[self.D11, self.A12], [self.A21, self.A22]])
            self.deltaQ_deltaH = -(np.linalg.inv(self.DA)) @ self.F
            # print('\nThe deltaQ deltaH Matrix in iteration: ' + str(count) + '\n')
            # print(self.deltaQ_deltaH)
            self.q_h += self.deltaQ_deltaH
            # print('\nThe updated QH Matrix in iteration: ' + str(count) + '\n')
            # print(self.q_h)
            self.D11 = np.diag(self.K * self.q_h[0:self.pipe_num].T[0] ** (self.n - 1))
            self.q_star = self.A21 @ self.q_h[0:self.pipe_num]
            self.b = np.concatenate([-self.A10 @ self.H0, self.q_star])
            self.F = np.concatenate(
                [self.D11 @ self.q_h[0:self.pipe_num] + self.A12 @ self.q_h[
                                                                   self.pipe_num:self.pipe_num + self.node_num],
                 self.A21 @ self.q_h[0:self.pipe_num]]) - self.b
            dE = -self.F[0:self.pipe_num]
            # print('\nThe updated dE Matrix in iteration: ' + str(count) + '\n')
            # print(dE)
            dq = -self.F[self.pipe_num:self.pipe_num + self.node_num]
            # print('\nThe updated dq Matrix in iteration: ' + str(count) + '\n')
            # print(dq)
            count += 1
        # print(count)
