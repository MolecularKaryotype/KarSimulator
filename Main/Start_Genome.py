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

                    full_KT[chr_name].append(Chromosome(chr_name, Arm([p_arm_segment], 'p'), Arm([q_arm_segment], 'q'),
                                                        t1_len, t2_len, Arm([centromere_segment], 'centromere')))

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
    initialization_string = 'initialization\n\tautosomes: {}\n\tautosomal copy number: {}\n\tsex chromosomes: {}' \
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
    ordinal_history = []

    # read-in segment indexing
    with open(input_file) as fp_read:
        line = fp_read.readline()  # column headers
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

        line = fp_read.readline()  # column headers

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
                new_chromosome = Chromosome(info[0], Arm([], 'deleted'), Arm([], 'deleted'), int(info[2]),
                                            int(info[3]), Arm([], 'deleted'), True)
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

            new_chromosome = Chromosome(info[0], Arm(p_arm_segments, 'p'), Arm(q_arm_segments, 'q'),
                                        int(info[2]), int(info[3]),
                                        Arm(current_centromere_segments, 'centromere'))
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
        block_name = 'null_block_error'
        first_block_passed = False
        while True:
            line = fp_read.readline().replace('\n', '')
            if not line:
                break
            if line[0] != '\t':
                # history block marking
                if first_block_passed:
                    history_markers[len(histories) - 1] = block_name
                block_name = line.split(': ')[1]
                first_block_passed = True
            else:
                line = line.split('\t')
                ordinal_history_info = line[2]
                line = line[1]
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

                # store current history entry
                history_segments = segment_list_to_segments(segment_list)
                history_chr_from = get_chromosome_from_unique_id(chr_from_name)
                history_chr_to = get_chromosome_from_unique_id(chr_to_name)
                histories.append(tuple([event_type, Arm(history_segments, 'history'),
                                        history_chr_from, history_chr_to]))

                # store current ordinal history entry
                ordinal_history_info = ordinal_history_info.split('(')
                if len(ordinal_history_info) == 1:
                    ordinal_history.append([])  # empty list if no ordinal data is given
                else:
                    new_ordinal_entry = []

                    for entry_index in range(1, len(ordinal_history_info)):
                        current_entry_info = ordinal_history_info[entry_index].replace(')', '').split(',')
                        current_entry = [current_entry_info[0]]
                        if current_entry_info[0] in ['inv', 'ins']:
                            # inv and ins are always chr_to
                            event_segment_index = current_entry_info[1].split('.')[0]
                            event_segment_ordinal = int(current_entry_info[1].split('.')[1])
                            event_segment = segment_dict[event_segment_index].duplicate()
                            event_segment.ordinal = event_segment_ordinal

                            is_on_q_arm, segment_arm_location = get_segment_location(event_segment, event_segment_ordinal, history_chr_to)
                            left_boundary_segment = None
                            right_boundary_segment = None
                            if not is_on_q_arm:
                                left_boundary_segment_location = segment_arm_location - 1
                                right_boundary_segment_location = segment_arm_location + 1
                                if left_boundary_segment_location >= 0:
                                    left_boundary_segment = history_chr_to.p_arm.segments[left_boundary_segment_location].duplicate()
                                if right_boundary_segment_location < len(history_chr_to.p_arm.segments):
                                    right_boundary_segment = history_chr_to.p_arm.segments[right_boundary_segment_location].duplicate()
                            else:
                                left_boundary_segment_location = segment_arm_location - 1
                                right_boundary_segment_location = segment_arm_location + 1
                                if left_boundary_segment_location >= 0:
                                    left_boundary_segment = history_chr_to.q_arm.segments[left_boundary_segment_location].duplicate()
                                if right_boundary_segment_location < len(history_chr_to.q_arm.segments):
                                    right_boundary_segment = history_chr_to.q_arm.segments[right_boundary_segment_location].duplicate()

                            if left_boundary_segment is not None:
                                left_boundary_segment.ordinal = get_segment_ordinal(left_boundary_segment, history_chr_to)
                            if right_boundary_segment is not None:
                                right_boundary_segment.ordinal = get_segment_ordinal(right_boundary_segment, history_chr_to)
                            current_entry += [event_segment, left_boundary_segment, right_boundary_segment]
                            new_ordinal_entry.append(tuple(current_entry))
                        elif current_entry_info[0] == 'del':
                            # deletion is always on chr_from
                            event_segment = history_segments
                            left_boundary_index = None
                            right_boundary_index = None
                            left_boundary_segment = current_entry_info[1]
                            right_boundary_segment = current_entry_info[2]
                            if current_entry_info[1] not in ['T1', 'T2', 'CEN']:
                                left_boundary_index = current_entry_info[1].split('.')[0]
                                left_boundary_ordinal = int(current_entry_info[1].split('.')[1])
                            if current_entry_info[2] not in ['T1', 'T2', 'CEN']:
                                right_boundary_index = current_entry_info[2].split('.')[0]
                                right_boundary_ordinal = int(current_entry_info[2].split('.')[1])

                            if left_boundary_index is not None:
                                left_boundary_segment = segment_dict[left_boundary_index].duplicate()
                                left_boundary_segment.ordinal = left_boundary_ordinal
                            if right_boundary_index is not None:
                                right_boundary_segment = segment_dict[right_boundary_index].duplicate()
                                right_boundary_segment.ordinal = right_boundary_ordinal
                            current_entry += [event_segment, left_boundary_segment, right_boundary_segment]
                            new_ordinal_entry.append(tuple(current_entry))

                    ordinal_history.append(new_ordinal_entry)

        if first_block_passed:
            history_markers[len(histories) - 1] = block_name  # append last block

    return Genome(full_KT, motherboard_segments, centromere_segments, initialization_string, histories, history_markers, ordinal_history=ordinal_history)


def get_segment_location(input_current_segment: Segment, input_segment_ordinal: int, input_chr: Chromosome):
    output_segment_location = -1
    finder_ptr = 0
    p_arm_exhausted = False
    for occurrence in range(0, input_segment_ordinal):
        segment_not_matched = True
        while segment_not_matched:
            if not p_arm_exhausted:
                if input_chr.p_arm.segments[finder_ptr].same_segment_ignore_dir(input_current_segment):
                    if occurrence + 1 == input_segment_ordinal:
                        output_segment_location = finder_ptr
                    else:
                        finder_ptr += 1
                        if finder_ptr >= len(input_chr.p_arm.segments):
                            p_arm_exhausted = True
                            finder_ptr = 0
                    segment_not_matched = False
                else:
                    finder_ptr += 1
                    if finder_ptr >= len(input_chr.p_arm.segments):
                        p_arm_exhausted = True
                        finder_ptr = 0
            else:
                if input_chr.q_arm.segments[finder_ptr].same_segment_ignore_dir(input_current_segment):
                    if occurrence + 1 == input_segment_ordinal:
                        output_segment_location = finder_ptr
                    else:
                        finder_ptr += 1
                        if finder_ptr >= len(input_chr.q_arm.segments):
                            raise RuntimeError('segment not found with the correct ordinal in given chr')
                    segment_not_matched = False
                else:
                    finder_ptr += 1
                    if finder_ptr >= len(input_chr.q_arm.segments):
                        raise RuntimeError('segment not found with the correct ordinal in given chr')
    return p_arm_exhausted, output_segment_location


def get_segment_ordinal(input_segment: Segment, input_chr: Chromosome):
    occurrence = 0
    for segment_itr in input_chr:
        # FIXME: never return True
        if input_segment is segment_itr:
            occurrence += 1
            break
        elif input_segment.same_segment_ignore_dir(segment_itr):
            occurrence += 1
    return occurrence


def test():
    genome = generate_genome_from_KT('/media/zhaoyang-new/workspace/KarSim/0922_test/r1/24XYe10_r1.kt.txt')
    genome.output_KT('/media/zhaoyang-new/workspace/KarSim/0922_test/r1/24XYe10_r1_copy.kt.txt')


def test_ordinal_update():
    genome = generate_genome_from_KT('/media/zhaoyang-new/workspace/KarSim/KarSimulator/test_folder/23Xe10_r1.kt.txt')
    print('x')


if __name__ == "__main__":
    test_ordinal_update()
