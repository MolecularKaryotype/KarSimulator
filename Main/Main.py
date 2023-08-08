from Structures import *


def generate_raw_genome(copy_number: int, autosomes: [str], sex_chromosomes: [str], genome_index_file: str) -> Genome:
    """
    for generating raw KT files containing only the chromosome of interest
    :param copy_number: integer >= 1
    :param autosomes: ALL for all chromosomes, else input ['Chr1', 'Chr12', etc.]
    :param sex_chromosomes: Male for ['ChrX', 'ChrY'], Female for ['ChrX', 'ChrX'], else input directly
    :param genome_index_file: file containing the positioning information of the genome,
    see Metadata/Full_Genome_Indices.txt for example
    :return: Genome Object
    """
    # construct KT slots
    full_KT = {}
    motherboard_segments = []
    centromere_segments = []

    # setup autosomes
    if autosomes[0].lower() == 'all':
        autosomes = ['Chr1', 'Chr2', 'Chr3', 'Chr4', 'Chr5', 'Chr6', 'Chr7', 'Chr8', 'Chr9', 'Chr10', 'Chr11',
                     'Chr12', 'Chr13', 'Chr14', 'Chr15', 'Chr16', 'Chr17', 'Chr18', 'Chr19', 'Chr20', 'Chr21',
                     'Chr22']
    for slot in autosomes:
        full_KT[slot] = []

    # setup sex chromosomes
    if sex_chromosomes[0].lower() == 'male':
        sex_chromosomes = ['ChrX', 'ChrY']
    elif sex_chromosomes[0].lower() == 'female':
        sex_chromosomes = ['ChrX', 'ChrX']
    for slot in set(sex_chromosomes):
        full_KT[slot] = []

    # prepare Raw KT
    def add_chromosome(chr_of_interest):
        with open(genome_index_file) as fp_read:
            for line in fp_read:
                line = line.replace('\n', '').split('\t')
                if line[0] == chr_of_interest:
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

    for autosome_chr_itr in autosomes:
        add_chromosome(autosome_chr_itr)
    for sex_chr_itr in sex_chromosomes:
        add_chromosome(sex_chr_itr)

    # setup copy numbers for autosomes
    if copy_number <= 0:
        raise ValueError('copy number')
    for slot in full_KT:
        if slot in autosomes:
            first_chromosome = full_KT[slot][0]
            for i in range(copy_number - 1):
                full_KT[slot].append(first_chromosome.duplicate())

    # rename each chromosome with unique name
    for slot in full_KT:
        current_slot = full_KT[slot]
        for i in range(len(current_slot)):
            current_slot[i].name = current_slot[i].name + chr(i + 97)

    # document initialization
    initialization_string = 'initialization\n\tautosomes: {}\n\tautosomal copy number: {}\n\tsex chromosomes: {}'\
        .format(str(autosomes), str(copy_number), str(sex_chromosomes))

    return Genome(full_KT, motherboard_segments, centromere_segments, initialization_string)


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

    # read-in segment indexing
    with open(input_file) as fp_read:
        line = fp_read.readline()   # column headers
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

        line = fp_read.readline()   # column headers

        # construct slots
        for chromosome in chromosome_of_interest:
            full_KT[chromosome] = []

        # read-in the segments in each chromosome
        while True:
            line = fp_read.readline().replace('\n', '')
            if line[0] == '-':
                break
            info = line.split('\t')
            p_arm_segments = []
            q_arm_segments = []
            current_centromere_segments = []
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
            current_centromere_segments.append(current_segment)
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
                                        Arm(current_centromere_segments))
            # append chromosome to slot in genome, ignore the last char in chromosome's unique ID
            full_KT[info[0][:-1]].append(new_chromosome)

        # read in initialization information
        line = fp_read.readline()  # initialization title
        initialization_string = line
        while True:
            line = fp_read.readline()
            if line[0] != '\t':
                break
            initialization_string += line

    return Genome(full_KT, motherboard_segments, centromere_segments, '')


# new_genome = generate_raw_genome(2, ['Chr1', 'Chr2'], ['ChrX', 'ChrX', 'ChrY'],
#                                  '../Metadata/test_Full_Genome_Indices.txt')
# chr_1a = new_genome.full_KT['Chr1'][0]
# chr_1b = new_genome.full_KT['Chr1'][1]
# chr_2a = new_genome.full_KT['Chr2'][0]
# chr_2b = new_genome.full_KT['Chr2'][1]
# new_genome.inversion(chr_1a, chr_1a.p_arm, 27, 53)
# new_genome.deletion(chr_1a, chr_1a.p_arm, 28, 29)
# new_genome.mark_history('test_1')
# new_genome.right_duplication_inversion(chr_2a, chr_2a.p_arm, 27, 53)
# new_genome.left_duplication_inversion(chr_2b, chr_2b.p_arm, 27, 53)
# new_genome.translocation_reciprocal(chr_1a, chr_1a.p_arm, 1, 20,
#                                     chr_2a, chr_2a.q_arm, 1, 20)
# new_genome.mark_history('test_2')
# new_genome.output_KT('RAW_test_genome4.txt')

# new_genome = generate_genome_from_KT('../RAW_KT_hg38.txt')
# new_genome = generate_genome_from_KT('./RAW_test_genome.txt')
# new_genome = generate_genome_from_KT('./RAW_test_genome4.txt')
# print(new_genome.motherboard_tostring())
# print(new_genome.history_tostring())
# print(new_genome.KT_tostring())

# new_genome = generate_raw_genome(2, ['ALL'], '../Metadata/Full_Genome_Indices.txt')
# new_genome.output_KT('RAW_full_genome.txt')
