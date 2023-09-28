from Main.Start_Genome import generate_raw_genome

new_genome = generate_raw_genome(2, ['all'], ['female'],
                                 '../Metadata/Full_Genome_Indices.txt')
chr_17a = new_genome.full_KT['Chr17'][0]
event_segments = new_genome.tandem_duplication(chr_17a.p_arm, 16809758, 20258836)
new_genome.append_history('tandem duplication', event_segments, chr_17a, chr_17a)

new_genome.mark_history('Potocki_Lupski')
new_genome.output_KT('../Precompiled_Kar/Potocki_Lupski.kt.txt')
