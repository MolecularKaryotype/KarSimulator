import sys
sys.path.insert(1, '/media/zhaoyang-new/workspace/KarSim/KarSimulator/Main')
from Main.Start_Genome import generate_raw_genome


new_genome = generate_raw_genome(1, ['all'], ['ChrX'],
                                 '../Metadata/Full_Genome_Indices.txt')
event_segments = new_genome.tandem_duplication(new_genome.full_KT['Chr1'][0].q_arm, 21877244, 23226635)
new_genome.append_history('tandem_duplication', event_segments, new_genome.full_KT['Chr1'][0], new_genome.full_KT['Chr1'][0])

new_genome.mark_history('1q21_recurrent_microduplication')
new_genome.output_KT('../Precompiled_Kar/23X_1q21_recurrent_microduplication.kt.txt')
