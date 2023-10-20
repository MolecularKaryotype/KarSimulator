from Main.Start_Genome import generate_raw_genome

new_genome = generate_raw_genome(1, ['all'], ['ChrX'],
                                 '../Metadata/Full_Genome_Indices.txt')
event_segments = new_genome.deletion(new_genome.full_KT['Chr22'][0].q_arm, 6508509, 8325939)
new_genome.append_history('deletion', event_segments, new_genome.full_KT['Chr22'][0], new_genome.full_KT['Chr22'][0])

new_genome.mark_history('22q11.2_distal_deletion')
new_genome.output_KT('../Precompiled_Kar/23X_22q11.2_distal_deletion.kt.txt')
