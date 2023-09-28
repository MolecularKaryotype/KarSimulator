from Main.Start_Genome import generate_raw_genome

new_genome = generate_raw_genome(2, ['all'], ['male'],
                                 '../Metadata/Full_Genome_Indices.txt')
chr_Ya = new_genome.full_KT['ChrY'][0]
event_segments = new_genome.deletion(chr_Ya.q_arm, 7462125, 13375010)
new_genome.append_history('deletion', event_segments, chr_Ya, chr_Ya)

new_genome.mark_history('AZFb')
new_genome.output_KT('../Precompiled_Kar/AZFb.kt.txt')
