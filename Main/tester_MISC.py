import IO
from Start_Genome import *

genome = generate_genome_from_KT('TDS2.txt')
chr_12a = genome.full_KT['Chr12'][0]

# for chromosome in genome:
#     print(chromosome.name)
#     for segment in chromosome:
#         print(segment)
print(genome)
print(genome.KT_tostring())