from Main.Start_Genome import generate_raw_genome

new_genome = generate_raw_genome(2, ['Chr1', 'Chr2'], ['male'],
                                 '../Metadata/Full_Genome_Indices.txt')

new_genome.output_KT('../demo_folder/demo_small.kt.txt')
