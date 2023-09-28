from Main.Start_Genome import generate_raw_genome

new_genome = generate_raw_genome(2, ['all'], ['female'],
                                 '../Metadata/Full_Genome_Indices.txt')
chr_11a = new_genome.full_KT['Chr11'][0]
event_segments = new_genome.deletion(chr_11a.p_arm, 43913250, 45970899)
new_genome.append_history('deletion', event_segments, chr_11a, chr_11a)

new_genome.mark_history('Potocki Shaffer')
new_genome.output_KT('../Precompiled_Kar/Potocki_Shaffer.kt.txt')
