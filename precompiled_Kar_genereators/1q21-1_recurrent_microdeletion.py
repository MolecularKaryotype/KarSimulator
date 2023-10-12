from Main.Start_Genome import generate_raw_genome

new_genome = generate_raw_genome(2, ['all'], ['female'],
                                 '../Metadata/Full_Genome_Indices.txt')
event_segments = new_genome.deletion(new_genome.full_KT['Chr1'][0].q_arm, 21877244, 23226635)
new_genome.append_history('deletion', event_segments, new_genome.full_KT['Chr1'][0], new_genome.full_KT['Chr1'][0])

new_genome.mark_history('1q21-1_recurrent_microdeletion')
new_genome.output_KT('../Precompiled_Kar/1q21-1_recurrent_microdeletion.kt.txt')
