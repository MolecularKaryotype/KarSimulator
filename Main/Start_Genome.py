from Main.Structures import *


def generate_raw_genome(copy_number: int, autosomes: [str], sex_chromosomes: [str], genome_index_file: str) -> Genome:
    """
    for generating raw KT files containing only the chromosome of interest
    :param copy_number: integer >= 1
    :param autosomes: ALL for all chromosomes, else input ['Chr1', 'Chr12', etc.]
    :param sex_chromosomes: Male for ['ChrX', 'ChrY'], Female for ['ChrX', 'ChrX'], else input directly
    :param genome_index_file: file containing the positioning information of the genome,
    see Metadata/hg38_index.txt for example
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
    if len(sex_chromosomes) > 0:
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

                    def segment_contained_in_list(item, this_list):
                        return any(item == e for e in this_list)

                    def add_if_no_duplicate(segment: Segment, this_list):
                        deep_copied_segment = segment.duplicate()
                        if deep_copied_segment not in this_list:
                            this_list.append(deep_copied_segment)

                    add_if_no_duplicate(p_arm_segment, motherboard_segments)
                    add_if_no_duplicate(q_arm_segment, motherboard_segments)
                    add_if_no_duplicate(centromere_segment, centromere_segments)

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
    histories = []
    history_markers = {}

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

            if info[1] == 'deleted':
                new_chromosome = Chromosome(info[0], Arm([]), Arm([]), int(info[2]),
                                            int(info[3]), Arm([]), True)
                full_KT[info[0][:-1]].append(new_chromosome)
                continue

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
            index_direction = segment_indices[current_index][-1]
            current_segment = segment_dict[segment_indices[current_index][:-1]].duplicate()
            if index_direction == "-":
                current_segment.invert()
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
        for i in range(3):
            line = fp_read.readline()
            if i == 2:
                initialization_string += line.replace('\n', '')
                break
            initialization_string += line

        # read in history
        block_name = 'null_block'
        first_block_passed = False
        while True:
            line = fp_read.readline().replace('\n', '')
            if not line:
                break
            if line[0] != '\t':
                if first_block_passed:
                    history_markers[len(histories) - 1] = block_name
                block_name = line.split(': ')[1]
                first_block_passed = True
            else:
                line = line.replace('\t', '')
                line = line.split('], ')
                event_type = line[0].split(' on')[0]
                segment_list = line[0].split('segments [')[1].split(',')
                line = line[1].replace('from ', '').split(' to ')
                chr_from_name = line[0]
                chr_to_name = line[1]

                def segment_list_to_segments(input_segment_list):
                    output_segments = []
                    for segment_index_reference in input_segment_list:
                        sign = segment_index_reference[-1]
                        segment_index_reference = segment_index_reference[:-1]
                        this_segment = segment_dict[segment_index_reference]
                        if sign == '+':
                            output_segments.append(this_segment)
                        else:
                            reversed_this_segment = this_segment.duplicate()
                            reversed_this_segment.invert()
                            output_segments.append(reversed_this_segment)
                    return output_segments

                def get_chromosome_from_unique_id(unique_id):
                    this_slot = unique_id[:-1]
                    current_slot = full_KT[this_slot]
                    for chromosome_itr in current_slot:
                        if chromosome_itr.name == unique_id:
                            return chromosome_itr
                    raise IndexError('get_chromosome_from_unique_id ID not found in Full_KT')

                history_segments = segment_list_to_segments(segment_list)
                history_chr_from = get_chromosome_from_unique_id(chr_from_name)
                history_chr_to = get_chromosome_from_unique_id(chr_to_name)
                histories.append(tuple([event_type, Arm(history_segments), history_chr_from, history_chr_to]))

        history_markers[len(histories) - 1] = block_name    # append last block

    return Genome(full_KT, motherboard_segments, centromere_segments, initialization_string, histories, history_markers)