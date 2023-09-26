from Main.Start_Genome import generate_raw_genome

new_genome = generate_raw_genome(2, ['all'], ['male'],
                                 '../Metadata/Full_Genome_Indices.txt')
chr_22a = new_genome.full_KT['Chr22'][0]
chr_22b = new_genome.full_KT['Chr22'][1]
event_segments = new_genome.tandem_duplication(chr_22a.p_arm, 0, 2444788)
new_genome.append_history('tandem duplication', event_segments, chr_22a, chr_22a)
event_segments = new_genome.tandem_duplication(chr_22b.p_arm, 0, 2444788)
new_genome.append_history('tandem duplication', event_segments, chr_22b, chr_22b)
event_segments = new_genome.tandem_duplication(chr_22a.q_arm, 0, 3473482)
new_genome.append_history('tandem duplication', event_segments, chr_22a, chr_22a)
event_segments = new_genome.tandem_duplication(chr_22b.q_arm, 0, 3473482)
new_genome.append_history('tandem duplication', event_segments, chr_22b, chr_22b)

new_genome.mark_history('Cat Eye 1')
new_genome.output_KT('../Precompiled_Kar/Cat_Eye_1.kt.txt')
