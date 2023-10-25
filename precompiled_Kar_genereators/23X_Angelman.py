import sys
sys.path.insert(1, '/media/zhaoyang-new/workspace/KarSim/KarSimulator/Main')
from Main.Start_Genome import generate_raw_genome


new_genome = generate_raw_genome(1, ['all'], ['ChrX'],
                                 '../Metadata/Full_Genome_Indices.txt')
event_segments = new_genome.deletion(new_genome.full_KT['Chr15'][0].q_arm, 2952090, 8467865)
new_genome.append_history('deletion', event_segments, new_genome.full_KT['Chr15'][0], new_genome.full_KT['Chr15'][0])

new_genome.mark_history('Angelman')
new_genome.output_KT('../Precompiled_Kar/23X_Angelman.kt.txt')
