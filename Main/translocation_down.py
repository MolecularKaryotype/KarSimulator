from Main import generate_raw_genome, generate_genome_from_KT

new_genome = generate_raw_genome(2, ['Chr12', 'Chr21'], ['male'],
                                 '../Metadata/Full_Genome_Indices.txt')
chr_12a = new_genome.full_KT['Chr12'][0]
chr_21a = new_genome.full_KT['Chr21'][0]

new_genome.chromosomal_duplication(chr_21a)
chr_21c = new_genome.full_KT['Chr21'][2]
new_genome.translocation_reciprocal(chr_12a, chr_12a.q_arm, 0, 96080055,
                                    chr_21c, chr_21c.q_arm, 0, 33784173)
new_genome.chromosomal_deletion(chr_21c)
new_genome.mark_history('Translocation Down Syndrome')
new_genome.output_KT('TDS1.txt')
# new_genome.output_FASTA('../Preparation/All24Chr.fasta', 'TDS1v4.fasta')

new_genome = generate_genome_from_KT('TDS1.txt')
chr_12a = new_genome.full_KT['Chr12'][0]
chr_12b = new_genome.full_KT['Chr12'][1]
chr_21a = new_genome.full_KT['Chr21'][0]
chr_Xa = new_genome.full_KT['ChrX'][0]
new_genome.deletion(chr_12a, chr_12a.p_arm, 14563405, 21563405)
new_genome.inversion(chr_12b, chr_12b.q_arm, 29563405, 73563405)
new_genome.duplication(chr_21a, chr_21a.q_arm, 12915809, 17915809)
new_genome.right_duplication_inversion(chr_Xa, chr_Xa.q_arm, 9915809, 25915809)
new_genome.mark_history('random mutations')
new_genome.output_KT('TDS2.txt')
new_genome.output_FASTA('../Preparation/All24Chr.fasta', 'TDS2v7.fasta')
