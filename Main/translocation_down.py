from Main import generate_raw_genome, generate_genome_from_KT

new_genome = generate_raw_genome(2, ['Chr14', 'Chr21'], ['male'],
                                 '../Metadata/Full_Genome_Indices.txt')
chr_14a = new_genome.full_KT['Chr14'][0]
chr_21a = new_genome.full_KT['Chr21'][0]

new_genome.chromosomal_duplication(chr_21a)
chr_21c = new_genome.full_KT['Chr21'][2]
new_genome.translocation_reciprocal(chr_14a, chr_14a.q_arm, 0, 88710193,
                                    chr_21c, chr_21c.q_arm, 0, 33784173)
new_genome.chromosomal_deletion(chr_21c)
new_genome.mark_history('Translocation Down Syndrome')
new_genome.output_KT('TDS1.txt')

