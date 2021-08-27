from GGA import GGA
from Partition import Partition
from Inputs import Inputs
import datetime

gga = GGA('./Benchmarks/1_Hanoi.inp')
print(gga.adjacency_list)

# start = datetime.datetime.now()
# gga.gga_algorithm()
# end = datetime.datetime.now()

# print('Time : ' + str((end - start).total_seconds() * 1000))

# partition = Partition()
# partition.partition()
# partition.get_connection()

# inputs = Inputs('./Benchmarks/1_Hanoi.inp')
# inputs.find_and_map()
# inputs.make_A()
