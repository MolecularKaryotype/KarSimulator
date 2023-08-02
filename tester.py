from Main.Structures import *

ALL_chr = ['Chr1', 'Chr2', 'Chr3', 'Chr4', 'Chr5', 'Chr6', 'Chr7', 'Chr8', 'Chr9', 'Chr10', 'Chr11', 'Chr12',
           'Chr13', 'Chr14', 'Chr15', 'Chr16', 'Chr17', 'Chr18', 'Chr19', 'Chr20', 'Chr21', 'Chr22', 'ChrX',
           'ChrY']

new_genome = Genome(2, ['Chr1', 'Chr2'], './Metadata/test_Full_Genome_Indices.txt')
chr_1a = new_genome.full_KT['Chr1'][0]
chr_1b = new_genome.full_KT['Chr1'][1]
chr_2a = new_genome.full_KT['Chr2'][0]
chr_2b = new_genome.full_KT['Chr2'][1]
chr_1a.p_arm.segments = [Segment('Chr1', 0, 25), Segment('Chr1', 26, 37),
                         Segment('Chr1', 38, 39), Segment('Chr1', 40, 50),
                         Segment('Chr1', 51, 76), Segment('Chr1', 77, 100)]
chr_1b.p_arm.segments = [Segment('Chr1', 0, 25), Segment('Chr1', 26, 37),
                         Segment('Chr1', 38, 39), Segment('Chr1', 40, 50),
                         Segment('Chr1', 51, 76), Segment('Chr1', 77, 100)]
# chr_2a.p_arm.segments = [Segment('Chr1', 0, 25), Segment('Chr1', 26, 37),
#                          Segment('Chr1', 38, 39), Segment('Chr1', 40, 50),
#                          Segment('Chr1', 51, 76), Segment('Chr1', 77, 100)]
# chr_2b.p_arm.segments = [Segment('Chr1', 0, 25), Segment('Chr1', 26, 37),
#                          Segment('Chr1', 38, 39), Segment('Chr1', 40, 50),
#                          Segment('Chr1', 51, 76), Segment('Chr1', 77, 100)]

new_genome.deletion(chr_1a, chr_1a.p_arm, 27, 53)
new_genome.deletion(chr_2a, chr_2a.p_arm, 27, 53)
print(new_genome)
