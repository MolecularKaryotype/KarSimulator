from Start_Genome import generate_raw_genome, generate_genome_from_KT
from Structures import *

# new_genome = generate_raw_genome(2, ['Chr1', 'Chr2'], ['ChrX', 'ChrX', 'ChrY'],
#                                  '../Metadata/test_Full_Genome_Indices.txt')
# chr_1a = new_genome.full_KT['Chr1'][0]
# chr_1b = new_genome.full_KT['Chr1'][1]
# chr_2a = new_genome.full_KT['Chr2'][0]
# chr_2b = new_genome.full_KT['Chr2'][1]
# new_genome.inversion(chr_1a, chr_1a.p_arm, 27, 53)
# new_genome.deletion(chr_1a, chr_1a.p_arm, 33, 35)
# new_genome.mark_history('manual_input')
# # new_genome.right_duplication_inversion(chr_2a, chr_2a.p_arm, 27, 53)
# # new_genome.left_duplication_inversion(chr_2b, chr_2b.p_arm, 27, 53)
# # new_genome.translocation_reciprocal(chr_1a, chr_1a.p_arm, 1, 20,
# #                                     chr_2a, chr_2a.q_arm, 1, 20)
# # new_genome.mark_history('test_2')
# new_genome.output_KT('test_read_KT_step1.txt')
#
# # new_genome = generate_genome_from_KT('../RAW_KT_hg38.txt')
# # new_genome = generate_genome_from_KT('./RAW_test_genome.txt')
# new_genome = generate_genome_from_KT('./test_read_KT_step1.txt')
# chr_1a = new_genome.full_KT['Chr1'][0]
# chr_1b = new_genome.full_KT['Chr1'][1]
# chr_2a = new_genome.full_KT['Chr2'][0]
# chr_2b = new_genome.full_KT['Chr2'][1]
# new_genome.right_duplication_inversion(chr_2a, chr_2a.p_arm, 27, 53)
# new_genome.left_duplication_inversion(chr_2b, chr_2b.p_arm, 27, 53)
# new_genome.translocation_reciprocal(chr_1a, chr_1a.p_arm, 1, 20,
#                                     chr_2a, chr_2a.q_arm, 1, 20)
# new_genome.mark_history('manual_input2')
# new_genome.output_KT('test_read_KT_step2.txt')
# # print(new_genome.motherboard_tostring())
# # print(new_genome.history_tostring())
# # print(new_genome.KT_tostring())
#
# # new_genome = generate_raw_genome(2, ['ALL'], '../Metadata/hg38_index.txt')
# # new_genome.output_KT('RAW_full_genome.txt')
#
# new_genome = generate_raw_genome(2, ['Chr1', 'Chr2'], ['ChrX', 'ChrX', 'ChrY'],
#                                  '../Metadata/test_Full_Genome_Indices.txt')
# chr_1a = new_genome.full_KT['Chr1'][0]
# chr_1b = new_genome.full_KT['Chr1'][1]
# chr_2a = new_genome.full_KT['Chr2'][0]
# chr_2b = new_genome.full_KT['Chr2'][1]
# new_genome.inversion(chr_1a, chr_1a.p_arm, 27, 53)
# new_genome.deletion(chr_1a, chr_1a.p_arm, 33, 35)
# new_genome.mark_history('manual_input')
# new_genome.right_duplication_inversion(chr_2a, chr_2a.p_arm, 27, 53)
# new_genome.left_duplication_inversion(chr_2b, chr_2b.p_arm, 27, 53)
# new_genome.translocation_reciprocal(chr_1a, chr_1a.p_arm, 1, 20,
#                                     chr_2a, chr_2a.q_arm, 1, 20)
# new_genome.mark_history('manual_input2')
# new_genome.output_KT('test_read_KT_together.txt')


## Testing chromosomal deletion
# new_genome = generate_raw_genome(2, ['Chr1', 'Chr2'], ['ChrX', 'ChrX', 'ChrY'],
#                                  '../Metadata/test_Full_Genome_Indices.txt')
# chr_1a = new_genome.full_KT['Chr1'][0]
# chr_1b = new_genome.full_KT['Chr1'][1]
# new_genome.inversion(chr_1a, chr_1a.p_arm, 27, 53)
# new_genome.mark_history('signature')
# new_genome.chromosomal_deletion(chr_1b)
# new_genome.mark_history('deletion1')
# new_genome.chromosomal_deletion(chr_1a)
# new_genome.mark_history('deletion2')
# new_genome.output_KT('test_chromosomal_deletion.txt')
#
# new_genome = generate_genome_from_KT('./test_chromosomal_deletion.txt')
# chr_2a = new_genome.full_KT['Chr2'][0]
# chr_2b = new_genome.full_KT['Chr2'][1]
# new_genome.chromosomal_deletion(chr_2a)
# new_genome.mark_history('deletion3')
# new_genome.output_KT('test_chromosomal_deletion2.txt')


## Testing chromosomal duplication
# new_genome = generate_raw_genome(2, ['Chr1', 'Chr2'], ['ChrX', 'ChrX', 'ChrY'],
#                                  '../Metadata/test_Full_Genome_Indices.txt')
# chr_1a = new_genome.full_KT['Chr1'][0]
# new_genome.chromosomal_duplication(chr_1a)
# chr_2a = new_genome.full_KT['Chr2'][0]
# new_genome.chromosomal_deletion(chr_2a)
# new_genome.mark_history('duplication')
# new_genome.output_KT('test_chromosomal_duplicate_delete.txt')
#
# new_genome = generate_genome_from_KT('./test_chromosomal_duplicate_delete.txt')
# new_genome.output_KT('test_chromosomal_duplicate_delete2.txt')


## Testing issue of translocation happening on the same arm
new_genome = generate_raw_genome(1, ['Chr1', 'Chr2', 'Chr3', 'Chr4'], ['ChrX', 'ChrY'],
                                 '../Metadata/Full_Genome_Indices.txt')
chr_2a = new_genome.full_KT['Chr2'][0]
new_genome.translocation_nonreciprocal(chr_2a.p_arm, 60000, 79999,
                                       chr_2a.p_arm, 20000)
new_genome.mark_history('trans')
new_genome.output_KT('/media/zhaoyang-new/workspace/KarSim/0922_test/test_nonreciprocal_trans.txt')
