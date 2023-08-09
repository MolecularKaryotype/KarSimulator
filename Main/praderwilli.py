from Main import generate_raw_genome, generate_genome_from_KT
from Structures import *

new_genome = generate_raw_genome(2, ['Chr1', 'Chr15'], ['ChrX', 'ChrX', 'ChrY'],
                                 '../Metadata/Full_Genome_Indices.txt')
chr_1a = new_genome.full_KT['Chr1'][0]
chr_15a = new_genome.full_KT['Chr15'][0]
chr_Xa = new_genome.full_KT['ChrX'][0]
new_genome.deletion(chr_15a, chr_15a.q_arm, 8000000, 8500000)
new_genome.mark_history('Prader Willi Signature')
new_genome.inversion(chr_15a, chr_15a.q_arm, 7000000, 9000000)
new_genome.deletion(chr_Xa, chr_Xa.q_arm, 1000000, 2000000)
new_genome.mark_history('Random Mutations1')
# new_genome.translocation_reciprocal(chr_15a, chr_15a.p_arm, 1000000, 2000000,
#                                     chr_15a, chr_15a.q_arm, 1000000, 2000000)
new_genome.translocation_reciprocal(chr_1a, chr_1a.p_arm, 100000, 200000,
                                    chr_15a, chr_15a.q_arm, 100000, 200000)
new_genome.mark_history('Random Mutations2')
new_genome.output_KT('PW_together.txt')

new_genome = generate_raw_genome(2, ['Chr1', 'Chr15'], ['ChrX', 'ChrX', 'ChrY'],
                                 '../Metadata/Full_Genome_Indices.txt')
chr_15a = new_genome.full_KT['Chr15'][0]
chr_Xa = new_genome.full_KT['ChrX'][0]
new_genome.deletion(chr_15a, chr_15a.q_arm, 8000000, 8500000)
new_genome.mark_history('Prader Willi Signature')
new_genome.output_KT('PW_step1.txt')
new_genome = generate_genome_from_KT('./PW_step1.txt')
chr_15a = new_genome.full_KT['Chr15'][0]
chr_Xa = new_genome.full_KT['ChrX'][0]
new_genome.inversion(chr_15a, chr_15a.q_arm, 7000000, 9000000)
new_genome.deletion(chr_Xa, chr_Xa.q_arm, 1000000, 2000000)
new_genome.mark_history('Random Mutations1')
new_genome.output_KT('PW_step2.txt')
new_genome = generate_genome_from_KT('./PW_step2.txt')
chr_1a = new_genome.full_KT['Chr1'][0]
chr_15a = new_genome.full_KT['Chr15'][0]
new_genome.translocation_reciprocal(chr_1a, chr_1a.p_arm, 100000, 200000,
                                    chr_15a, chr_15a.q_arm, 100000, 200000)
new_genome.mark_history('Random Mutations2')
new_genome.output_KT('PW_step3.txt')
