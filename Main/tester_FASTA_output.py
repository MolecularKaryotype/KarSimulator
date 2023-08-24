from Start_Genome import generate_raw_genome, generate_genome_from_KT
from Testers import test_KT_FASTA_Correspondence

# new_genome = generate_raw_genome(2, ['Chr12', 'Chr21'], [],
#                                  '../Metadata/test_Full_Genome_Indices.txt')
#
# chr_12a = new_genome.full_KT['Chr12'][0]
# chr_21a = new_genome.full_KT['Chr21'][0]
#
# new_genome.chromosomal_duplication(chr_21a)
# chr_21c = new_genome.full_KT['Chr21'][2]
# new_genome.translocation_reciprocal(chr_12a, chr_12a.q_arm, 0, 79,
#                                     chr_21c, chr_21c.q_arm, 0, 79)
# new_genome.chromosomal_deletion(chr_21c)
# new_genome.mark_history('Translocation Down Syndrome')
#
# chr_12a = new_genome.full_KT['Chr12'][0]
# chr_12b = new_genome.full_KT['Chr12'][1]
# chr_21a = new_genome.full_KT['Chr21'][0]
# chr_21b = new_genome.full_KT['Chr21'][1]
# new_genome.deletion(chr_12a, chr_12a.p_arm, 20, 30)
# new_genome.inversion(chr_12b, chr_12b.q_arm, 30, 40)
# new_genome.duplication(chr_21a, chr_21a.q_arm, 20, 30)
# new_genome.right_duplication_inversion(chr_21b, chr_21b.q_arm, 20, 40)
# new_genome.mark_history('random mutations')
#
# new_genome.output_KT('test_KT.txt')
# new_genome.output_FASTA('../Genomes/test_genome.fasta', 'test_fasta.fasta')

# genome = generate_genome_from_KT('test_KT.txt')
# test_KT_FASTA_Correspondence('test_KT.txt', 'test_fasta.fasta', '../Genomes/test_genome.fasta')



# new_genome = generate_raw_genome(2, ['Chr1'], [],
#                                  '../Metadata/test_Full_Genome_Indices2.txt')
# chr_1a = new_genome.full_KT['Chr1'][0]
# chr_1b = new_genome.full_KT['Chr1'][1]
# # new_genome.inversion(chr_1a, chr_1a.p_arm, 5, 12)
# # new_genome.inversion(chr_1b, chr_1b.q_arm, 0, 7)
# new_genome.inversion(chr_1a, chr_1a.p_arm, 1, 3)
# new_genome.inversion(chr_1b, chr_1b.q_arm, 1, 3)
# new_genome.mark_history('random mutations')
#
# new_genome.output_KT('test_KT2.txt')
# new_genome.output_FASTA('../Genomes/test_genome2.fasta', 'test_fasta2.fasta')

genome = generate_genome_from_KT('test_KT2.txt')
test_KT_FASTA_Correspondence('test_KT2.txt', 'test_fasta2.fasta', '../Genomes/test_genome2.fasta')
