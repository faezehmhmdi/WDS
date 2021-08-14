import numpy as np
import pymetis


class Partition:
    def __init__(self):
        self.num_of_subsystems = 2
        self.adjacency_list = [np.array([1, 2]),
                               np.array([0, 2]),
                               np.array([0, 1, 3]),
                               np.array([2, 4, 5]),
                               np.array([3, 5]),
                               np.array([3, 4])]
        self.n_cuts = 0
        self.membership = []

    def partition(self):
        self.n_cuts, self.membership = pymetis.part_graph(self.num_of_subsystems, adjacency=self.adjacency_list)
        # print(self.n_cuts)
        print(self.membership)
        for i in range(self.num_of_subsystems):
            print(np.argwhere(np.array(self.membership) == i).ravel())

    def get_connection(self):
        connections = []
        for i in range(self.num_of_subsystems):
            for j in np.argwhere(np.array(self.membership) == i).ravel():
                for k in self.adjacency_list[j]:
                    if self.membership[k] > self.membership[j]:
                        connections.append([k, j])
        print(connections)

    def make_matrix(self):
        pass
