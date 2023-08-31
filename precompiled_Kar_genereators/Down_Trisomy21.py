from Main.Start_Genome import generate_raw_genome

new_genome = generate_raw_genome(2, ['all'], ['male'],
                                 '../Metadata/Full_Genome_Indices.txt')
chr_21a = new_genome.full_KT['Chr21'][0]
new_genome.chromosomal_duplication(chr_21a)

new_genome.mark_history('Trisomy 21 Down Syndrome')
new_genome.output_KT('../Precompiled_Kar/Down_Trisomy21.kt.txt')
