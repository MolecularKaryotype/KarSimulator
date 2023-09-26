from Main.Start_Genome import generate_raw_genome

new_genome = generate_raw_genome(2, ['all'], ['male'],
                                 '../Metadata/Full_Genome_Indices.txt')
chr_5a = new_genome.full_KT['Chr5'][0]
event_segments = new_genome.deletion(chr_5a.p_arm, 0, 12523192)
new_genome.append_history('deletion', event_segments, chr_5a, chr_5a)

new_genome.mark_history('Cri du Chat')
new_genome.output_KT('../Precompiled_Kar/Cri_du_Chat.kt.txt')
