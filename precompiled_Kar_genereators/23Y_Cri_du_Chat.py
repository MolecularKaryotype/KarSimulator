from Main.Start_Genome import generate_raw_genome

new_genome = generate_raw_genome(1, ['all'], ['ChrY'],
                                 '../Metadata/Full_Genome_Indices.txt')
event_segments = new_genome.deletion(new_genome.full_KT['Chr5'][0].p_arm, 1, 12523192)
new_genome.append_history('deletion', event_segments, new_genome.full_KT['Chr5'][0], new_genome.full_KT['Chr5'][0])

new_genome.mark_history('Cri_du_Chat')
new_genome.output_KT('../Precompiled_Kar/23Y_Cri_du_Chat.kt.txt')
