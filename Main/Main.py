from Structures import *


def generate_raw_genome(copy_number: int, chr_of_interest: [str], genome_index_file: str) -> Genome:
    """
    for generating raw KT files containing only the chromosome of interest
    :param copy_number: integer >= 1
    :param chr_of_interest: ALL for all chromosomes, else input ['Chr1', 'ChrX', etc.]
    :param genome_index_file: file containing the positioning information of the genome,
    see Metadata/Full_Genome_Indices.txt for example
    :return: Genome Object
    """
    # construct KT slots
    full_KT = {}
    motherboard_segments = []
    centromere_segments = []

    if chr_of_interest[0] == 'ALL':
        chr_of_interest = ['Chr1', 'Chr2', 'Chr3', 'Chr4', 'Chr5', 'Chr6', 'Chr7', 'Chr8', 'Chr9', 'Chr10', 'Chr11',
                           'Chr12', 'Chr13', 'Chr14', 'Chr15', 'Chr16', 'Chr17', 'Chr18', 'Chr19', 'Chr20', 'Chr21',
                           'Chr22', 'ChrX', 'ChrY']
    for slot in chr_of_interest:
        full_KT[slot] = []
    # prepare Raw KT
    with open(genome_index_file) as fp_read:
        for line in fp_read:
            line = line.replace('\n', '').split('\t')
            if line[0] in chr_of_interest:
                chr_name = line[0]
                p_arm_segment = Segment(chr_name, int(line[2]), int(line[3]))
                q_arm_segment = Segment(chr_name, int(line[4]), int(line[5]))
                t1_len = int(line[2])
                t2_len = int(line[1]) - int(line[5]) - 1
                centromere_segment = Segment(chr_name, int(line[3]) + 1, int(line[4]) - 1)

                full_KT[chr_name].append(Chromosome(chr_name, Arm([p_arm_segment]), Arm([q_arm_segment]),
                                                    t1_len, t2_len, Arm([centromere_segment])))
                motherboard_segments.extend([p_arm_segment.duplicate(),
                                             q_arm_segment.duplicate()])
                centromere_segments.append(centromere_segment.duplicate())

    # setup copy numbers
    if copy_number <= 0:
        raise ValueError('copy number')
    for slot in full_KT:
        first_chromosome = full_KT[slot][0]
        for i in range(copy_number - 1):
            full_KT[slot].append(first_chromosome.duplicate())
    # rename each chromosome with unique name
    for slot in full_KT:
        current_slot = full_KT[slot]
        for i in range(copy_number):
            current_slot[i].name = current_slot[i].name + chr(i + 97)

    return Genome(full_KT, motherboard_segments, centromere_segments)


def generate_genome_from_KT(input_file: str) -> Genome:
    """
    directly generating Genome from a formatted KT file
    :param input_file: path to input file
    :return: Genome Object
    """
    full_KT = {}
    motherboard_segments = []
    centromere_segments = []
    chromosome_of_interest = set()
    segment_dict = {}

    with open(input_file) as fp_read:
        line = fp_read.readline()
        while True:
            line = fp_read.readline()
            if line[0] == '-':
                break
            info = line.replace('\n', '').split('\t')
            new_segment = Segment(info[1], int(info[2]), int(info[3]))
            segment_dict[info[0]] = new_segment
            chromosome_of_interest.add(info[1])
            if info[0][0] == "C":
                centromere_segments.append(new_segment)
            else:
                motherboard_segments.append(new_segment)

        line = fp_read.readline()
        # construct slots
        for chromosome in chromosome_of_interest:
            full_KT[chromosome] = []
        while True:
            line = fp_read.readline().replace('\n', '')
            if not line:
                break
            info = line.split('\t')
            p_arm_segments = []
            q_arm_segments = []
            centromere_segments = []
            segment_indices = info[1].split(',')
            # p_arm
            current_index = 0
            while segment_indices[current_index][0] != 'C':
                index_value = int(segment_indices[current_index][:-1])
                index_direction = segment_indices[current_index][-1]
                current_segment = segment_dict[str(index_value)].duplicate()
                if index_direction == "-":
                    current_segment.invert()
                p_arm_segments.append(current_segment)
                current_index += 1
            # centromere
            current_segment = segment_dict[segment_indices[current_index]].duplicate()
            centromere_segments.append(current_segment)
            current_index += 1
            # q_arm
            while current_index < len(segment_indices):
                index_value = int(segment_indices[current_index][:-1])
                index_direction = segment_indices[current_index][-1]
                current_segment = segment_dict[str(index_value)].duplicate()
                if index_direction == "-":
                    current_segment.invert()
                q_arm_segments.append(current_segment)
                current_index += 1

            new_chromosome = Chromosome(info[0], Arm(p_arm_segments), Arm(q_arm_segments), int(info[2]), int(info[3]),
                                        Arm(centromere_segments))
            full_KT[info[0][:-1]].append(new_chromosome)

    return Genome(full_KT, motherboard_segments, centromere_segments)


new_genome = generate_raw_genome(2, ['Chr1', 'Chr2'], '../Metadata/test_Full_Genome_Indices.txt')
chr_1a = new_genome.full_KT['Chr1'][0]
chr_1b = new_genome.full_KT['Chr1'][1]
chr_2a = new_genome.full_KT['Chr2'][0]
chr_2b = new_genome.full_KT['Chr2'][1]
new_genome.inversion(chr_1a, chr_1a.p_arm, 27, 53)
new_genome.deletion(chr_1a, chr_1a.p_arm, 28, 29)
new_genome.right_duplication_inversion(chr_2a, chr_2a.p_arm, 27, 53)
new_genome.left_duplication_inversion(chr_2b, chr_2b.p_arm, 27, 53)
new_genome.translocation_reciprocal(chr_1a, chr_1a.p_arm, 1, 20,
                                    chr_2a, chr_2a.q_arm, 1, 20)
new_genome.output_KT('RAW_test_genome.txt')

# new_genome = generate_genome_from_KT('../RAW_KT_hg38.txt')
# print(new_genome.motherboard_tostring())
# print(new_genome.history_tostring())
# print(new_genome.KT_tostring())

new_genome = generate_raw_genome(2, ['ALL'], '../Metadata/Full_Genome_Indices.txt')
new_genome.output_KT('RAW_full_genome.txt')
