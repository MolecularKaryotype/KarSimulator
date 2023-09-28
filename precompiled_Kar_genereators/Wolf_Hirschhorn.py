from Main.Start_Genome import generate_raw_genome

new_genome = generate_raw_genome(2, ['all'], ['male'],
                                 '../Metadata/Full_Genome_Indices.txt')
chr_4a = new_genome.full_KT['Chr4'][0]
event_segments = new_genome.deletion(chr_4a.p_arm, 1566470, 2098509)
new_genome.append_history('deletion', event_segments, chr_4a, chr_4a)


new_genome.mark_history('Wolf Hirschhorn')
new_genome.output_KT('../Precompiled_Kar/Wolf_Hirschhorn.kt.txt')
