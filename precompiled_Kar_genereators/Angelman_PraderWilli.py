from Main.Start_Genome import generate_raw_genome

# TODO: make a hg38 converter to convert hg38 index to p/q-arm indexing, using a metadata sheet as input (maybe set as
#  global variable)
new_genome = generate_raw_genome(2, ['all'], ['male'],
                                 '../Metadata/Full_Genome_Indices.txt')
chr_15a = new_genome.full_KT['Chr15'][0]
event_segments = new_genome.deletion(chr_15a.q_arm, 2952090, 8467865)
new_genome.append_history('deletion', event_segments, chr_15a, chr_15a)

new_genome.mark_history('Angelman PraderWilli')
new_genome.output_KT('../Precompiled_Kar/Angelman_PraderWilli.kt.txt')
