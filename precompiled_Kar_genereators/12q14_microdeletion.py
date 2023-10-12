from Main.Start_Genome import generate_raw_genome

new_genome = generate_raw_genome(2, ['all'], ['male'],
                                 '../Metadata/Full_Genome_Indices.txt')
event_segments = new_genome.deletion(new_genome.full_KT['Chr12'][0].q_arm, 27492886, 31066492)
new_genome.append_history('deletion', event_segments, new_genome.full_KT['Chr12'][0], new_genome.full_KT['Chr12'][0])

new_genome.mark_history('12q14_microdeletion')
new_genome.output_KT('../Precompiled_Kar/12q14_microdeletion.kt.txt')
