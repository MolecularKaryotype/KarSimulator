from Main.Start_Genome import generate_raw_genome

new_genome = generate_raw_genome(1, ['all'], ['ChrX'],
                                 '../Metadata/Full_Genome_Indices.txt')

new_genome.output_KT('../Precompiled_Kar/23X_ref.kt.txt')
