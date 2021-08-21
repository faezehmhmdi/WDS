from GGA import GGA
from Partition import Partition
from Inputs import Inputs

gga = GGA('./Benchmarks/1_Hanoi.inp')
# print(gga.A)
gga.gga_algorithm()
# partition = Partition()
# partition.partition()
# partition.get_connection()

# inputs = Inputs('./Benchmarks/1_Hanoi.inp')
# inputs.find_and_map()
# inputs.make_A()
