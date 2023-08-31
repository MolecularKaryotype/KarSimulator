from Main.Start_Genome import generate_raw_genome

new_genome = generate_raw_genome(2, ['all'], ['male'],
                                 '../Metadata/Full_Genome_Indices.txt')
chr_Xa = new_genome.full_KT['ChrX'][0]
new_genome.chromosomal_duplication(chr_Xa)

new_genome.mark_history('Klinefelter Syndrome')
new_genome.output_KT('../Precompiled_Kar/Klinefelter.kt.txt')
