from Main.Start_Genome import generate_raw_genome

new_genome = generate_raw_genome(2, ['all'], ['male'],
                                 '../Metadata/Full_Genome_Indices.txt')
chr_22a = new_genome.full_KT['Chr22'][0]
event_segments = new_genome.chromosomal_duplication(chr_22a)
chr_22c = new_genome.full_KT['Chr22'][2]
new_genome.append_history("chromosomal duplication", event_segments, chr_22a, chr_22c)
event_segments = new_genome.deletion(chr_22c.q_arm, 3473482, 35754148)
new_genome.append_history("deletion", event_segments, chr_22c, chr_22c)
event_segments = new_genome.chromosomal_duplication(chr_22c)
chr_22d = new_genome.full_KT['Chr22'][3]
new_genome.append_history("chromosomal duplication", event_segments, chr_22c, chr_22d)

new_genome.mark_history('Cat Eye 2')
new_genome.output_KT('../Precompiled_Kar/Cat_Eye_2.kt.txt')
