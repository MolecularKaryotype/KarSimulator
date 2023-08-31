from Main.Start_Genome import generate_raw_genome

new_genome = generate_raw_genome(2, ['all'], ['female'],
                                 '../Metadata/Full_Genome_Indices.txt')
new_genome.output_KT('../Precompiled_Kar/Female.kt.txt')
