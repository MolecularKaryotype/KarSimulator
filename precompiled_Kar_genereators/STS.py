from Main.Start_Genome import generate_raw_genome

new_genome = generate_raw_genome(2, ['all'], ['male'],
                                 '../Metadata/Full_Genome_Indices.txt')
chr_Xa = new_genome.full_KT['ChrX'][0]
event_segments = new_genome.deletion(chr_Xa.p_arm, 6527771, 8155154)
new_genome.append_history('deletion', event_segments, chr_Xa, chr_Xa)

new_genome.mark_history('STS')
new_genome.output_KT('../Precompiled_Kar/STS.kt.txt')
