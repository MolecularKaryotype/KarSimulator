import sys
sys.path.insert(1, '/media/zhaoyang-new/workspace/KarSim/KarSimulator/Main')
from Main.Start_Genome import generate_raw_genome


new_genome = generate_raw_genome(1, ['all'], ['ChrY'],
                                 '../Metadata/Full_Genome_Indices.txt')
event_segments = new_genome.tandem_duplication(new_genome.full_KT['Chr17'][0].p_arm, 14134598, 15507589)
new_genome.append_history('duplication', event_segments, new_genome.full_KT['Chr17'][0], new_genome.full_KT['Chr17'][0])

new_genome.mark_history('CMT1A')
new_genome.output_KT('../Precompiled_Kar/23Y_CMT1A.kt.txt')
