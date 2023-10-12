from Main.Start_Genome import generate_raw_genome

new_genome = generate_raw_genome(2, ['all'], ['male'],
                                 '../Metadata/Full_Genome_Indices.txt')
event_segments = new_genome.tandem_duplication(new_genome.full_KT['Chr17'][0].p_arm, 14194598, 15567589)
new_genome.append_history('duplication', event_segments, new_genome.full_KT['Chr17'][0], new_genome.full_KT['Chr17'][0])

new_genome.mark_history('CMT1A')
new_genome.output_KT('../Precompiled_Kar/CMT1A.kt.txt')
