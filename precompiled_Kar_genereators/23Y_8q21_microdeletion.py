import sys
sys.path.insert(1, '/media/zhaoyang-new/workspace/KarSim/KarSimulator/Main')
from Main.Start_Genome import generate_raw_genome


new_genome = generate_raw_genome(1, ['all'], ['ChrY'],
                                 '../Metadata/Full_Genome_Indices.txt')
event_segments = new_genome.deletion(new_genome.full_KT['Chr8'][0].q_arm, 30436963, 30976737)
new_genome.append_history('deletion', event_segments, new_genome.full_KT['Chr8'][0], new_genome.full_KT['Chr8'][0])

new_genome.mark_history('8q21_recurrent_microdeletion')
new_genome.output_KT('../Precompiled_Kar/23Y_8q21_recurrent_microdeletion.kt.txt')
