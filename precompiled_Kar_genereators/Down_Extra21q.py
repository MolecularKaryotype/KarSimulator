from Main.Start_Genome import generate_raw_genome

new_genome = generate_raw_genome(2, ['all'], ['male'],
                                 '../Metadata/Full_Genome_Indices.txt')
chr_12a = new_genome.full_KT['Chr12'][0]
chr_21a = new_genome.full_KT['Chr21'][0]

new_genome.chromosomal_duplication(chr_21a)
chr_21c = new_genome.full_KT['Chr21'][2]
new_genome.translocation_reciprocal_balanced(chr_12a.q_arm, 0, 96080055,
                                             chr_21c.q_arm, 0, 33784173)
new_genome.chromosomal_deletion(chr_21c)
new_genome.mark_history('Translocation Down Syndrome')
new_genome.output_KT('Down_Extra21q.kt.txt')
