from Main.Start_Genome import generate_raw_genome

new_genome = generate_raw_genome(1, ['all'], ['ChrX'],
                                 '../Metadata/Full_Genome_Indices.txt')
event_segments = new_genome.tandem_duplication(new_genome.full_KT['Chr15'][0].q_arm, 79089486, 82255934)
new_genome.append_history('tandem_duplication', event_segments, new_genome.full_KT['Chr15'][0], new_genome.full_KT['Chr15'][0])

new_genome.mark_history('15q26_overgrowth')
new_genome.output_KT('../Precompiled_Kar/23X_15q26_overgrowth.kt.txt')
