from Main.Start_Genome import generate_raw_genome

new_genome = generate_raw_genome(1, ['all'], ['ChrY'],
                                 '../Metadata/Full_Genome_Indices.txt')
event_segments = new_genome.deletion(new_genome.full_KT['Chr2'][0].p_arm, 59048561, 61582680)
new_genome.append_history('deletion', event_segments, new_genome.full_KT['Chr2'][0], new_genome.full_KT['Chr2'][0])

new_genome.mark_history('2p15-16.1_microdeletion')
new_genome.output_KT('../Precompiled_Kar/23Y_2p15-16.1_microdeletion.kt.txt')
