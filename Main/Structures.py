import IO


class Segment:
    chr_name: str
    start: int
    end: int

    def __init__(self, chr_name: str, start: int, end: int):
        self.chr_name = chr_name
        self.start = start
        self.end = end

    def __len__(self):
        return abs(self.end - self.start) + 1

    def __lt__(self, other):
        def get_chr_order(chromosome_name):
            chr_extracted = chromosome_name.replace('Chr', '')
            if chr_extracted == 'X':
                return 23
            elif chr_extracted == 'Y':
                return 24
            else:
                return int(chr_extracted)

        if get_chr_order(self.chr_name) < get_chr_order(other.chr_name):
            return True
        elif get_chr_order(self.chr_name) > get_chr_order(other.chr_name):
            return False
        elif get_chr_order(self.chr_name) == get_chr_order(other.chr_name):
            return max(self.start, self.end) < max(other.start, other.end)

    def __eq__(self, other):
        if isinstance(other, Segment):
            return (self.chr_name, self.start, self.end) == (other.chr_name, other.start, other.end)
        return False

    def __hash__(self):
        return hash((self.chr_name, self.start, self.end))

    def __str__(self):
        return "({}, {}, {})".format(self.chr_name, self.start, self.end)

    def same_segment_ignore_dir(self, other):
        if self.chr_name != other.chr_name:
            return False
        if self.start != other.start and self.start != other.end:
            return False
        if self.end != other.start and self.end != other.end:
            return False
        return True

    def to_string_ignore_dir(self):
        if self.direction():
            return "{}:{}-{}".format(self.chr_name, self.start, self.end)
        else:
            return "{}:{}-{}".format(self.chr_name, self.end, self.start)

    def direction(self):
        """
        :return: 1 for +, 0 for -
        """
        return self.start <= self.end

    def duplicate(self):
        return Segment(self.chr_name, self.start, self.end)

    def left_delete(self, bp_to_delete):
        """
        :param bp_to_delete: number of bp deleting
        :return: None
        """
        if self.direction():
            self.start = self.start + bp_to_delete
        else:
            self.start = self.start - bp_to_delete

    def right_delete(self, bp_to_delete):
        """
        :param bp_to_delete: number of bp deleting
        :return: None
        """
        if self.direction():
            self.end = self.end - bp_to_delete
        else:
            self.end = self.end + bp_to_delete

    def invert(self):
        temp_start = self.start
        self.start = self.end
        self.end = temp_start


class Arm:
    segments: [Segment]

    def __init__(self, segments: [Segment]):
        self.segments = segments

    def __len__(self):
        current_sum = 0
        for segment in self.segments:
            current_sum += len(segment)
        return current_sum

    def __contains__(self, item):
        return any(item == e for e in self.segments)

    def __str__(self):
        return_str = ''
        for segment in self.segments:
            return_str += str(segment)
        return return_str

    def duplicate(self):
        new_segments = []
        for segment in self.segments:
            new_segments.append(segment.duplicate())
        return Arm(new_segments)

    def delete_segments_by_index(self, segment_indices):
        self.segments = [segment for index, segment in enumerate(self.segments) if index not in segment_indices]

    def duplicate_segments_by_index(self, segment_indices):
        # segments come in order, so insert before the first segment
        index_of_insertion = segment_indices[0]
        new_segments = []
        for index in segment_indices:
            new_segments.append(self.segments[index].duplicate())

        self.segments[index_of_insertion:index_of_insertion] = new_segments

    def invert_segments_by_index(self, segment_indices):
        # segments come in order
        index_of_insertion = segment_indices[0]
        new_segments = []
        for index in reversed(segment_indices):
            new_segment = self.segments[index].duplicate()
            new_segment.invert()
            new_segments.append(new_segment)
        self.delete_segments_by_index(segment_indices)
        self.segments[index_of_insertion:index_of_insertion] = new_segments


class Chromosome:
    name: str
    p_arm: Arm
    q_arm: Arm
    centromere: Arm
    t1_len: int
    t2_len: int
    deleted: bool

    def __init__(self, name: str, p_arm: Arm, q_arm: Arm, t1_len: int, t2_len: int, centromere: Arm, deleted=False):
        self.name = name
        self.p_arm = p_arm
        self.q_arm = q_arm
        self.t1_len = t1_len
        self.t2_len = t2_len
        self.centromere = centromere
        self.deleted = deleted

    def __len__(self):
        return self.p_arm_len() + self.q_arm_len()

    def __str__(self):
        if self.deleted:
            return '{}: deleted'.format(self.name)
        return_str = '{}: t1-{} t2-{}\n\tp-arm: {}\n\tq-arm: {}\n\tCEN: {}' \
            .format(self.name, self.t1_len, self.t2_len, str(self.p_arm), str(self.q_arm), str(self.centromere))
        return return_str

    def __iter__(self):
        class ChromosomeIterator:
            def __init__(self, chromosome: Chromosome):
                self.chromosome = chromosome
                self.arms = [chromosome.p_arm, chromosome.centromere, chromosome.q_arm]
                self.current_arm_index = 0
                self.current_segment_index = 0

            def __next__(self):
                if self.current_arm_index < len(self.arms):
                    current_arm = self.arms[self.current_arm_index]
                    if self.current_segment_index < len(current_arm.segments):
                        segment = current_arm.segments[self.current_segment_index]
                        self.current_segment_index += 1
                        return segment
                    else:
                        self.current_arm_index += 1
                        self.current_segment_index = 0
                        return next(self)
                else:
                    raise StopIteration
        return ChromosomeIterator(self)

    def get_segment_from_range(self, left_bound, right_bound):
        """
        return a list of segments that are between (inclusive) the two bounds on this Chromosome object
        :param left_bound:
        :param right_bound:
        :return:
        """
        pass

    def p_arm_len(self):
        return len(self.p_arm)

    def q_arm_len(self):
        return len(self.q_arm)

    def duplicate(self):
        return Chromosome(self.name, self.p_arm.duplicate(), self.q_arm.duplicate(),
                          self.t1_len, self.t2_len, self.centromere.duplicate())


def index_global_to_arm(chromosome: Chromosome, left_index: int, right_index: int) -> (Chromosome, Arm, int, int):
    """
    given an absolute index range on a chromosome, locate the arm the range is on, and return the relative index range
    on that arm
    the index range cannot span different sections of a chromosome (sections: t1, p-arm, centromere, q-arm, t2)
    currently only support the index range that is on either p-arm or q-arm
    :param chromosome:
    :param left_index:
    :param right_index:
    :return: chromosome, arm on that chromosome, relative-left-index, relative-right-index
    """
    centromere_left_bound = chromosome.t1_len + chromosome.p_arm_len()
    centromere_right_bound = chromosome.t1_len + chromosome.p_arm_len() + len(chromosome.centromere) - 1
    if right_index < left_index:
        raise ValueError('absolute-left-index greater than absolute-right_index')
    elif left_index < chromosome.t1_len:
        raise ValueError('absolute-left-index in t1 region')
    elif right_index > len(chromosome) + len(chromosome.centromere):
        raise ValueError('absolute-right-index in t2 region')
    elif centromere_left_bound <= left_index <= centromere_right_bound:
        raise ValueError('absolute-left-index in centromere region')
    elif centromere_left_bound <= right_index <= centromere_right_bound:
        raise ValueError('absolute-right-index in centromere region')

    # locate the Arm left_index is on
    pass


class Genome:
    full_KT: {str: [Chromosome]}  # has as many slots as there are chromosome type, i.e. 24 for a male, 23 for a female
    motherboard: Arm  # using the Arm object to use generate breakpoint method
    centromere_segments = [Segment]
    history_block_markings = {}  # history enumerate index: block name
    history: [(str, Arm, Chromosome, Chromosome)]  # event type, event segments, chr from, chr to
    initialization_string: str  # contains information with the initialization of the genome

    def __init__(self, full_KT, motherboard_segments, centromere_segments, initialization_string,
                 history=None, history_markers=None):
        self.full_KT = full_KT
        self.motherboard = Arm(motherboard_segments)
        self.centromere_segments = centromere_segments
        self.initialization_string = initialization_string
        if history is not None:
            self.history = history
        else:
            self.history = []
        if history_markers is not None:
            self.history_block_markings = history_markers
        else:
            self.history_block_markings = {}

    def __str__(self):
        return_str = ''
        for chromosome in self:
            return_str += str(chromosome) + '\n'
        return return_str

    def __iter__(self):
        def custom_sort_chr(key):
            chr_part = key[3:]  # Extract the part after "Chr"
            if chr_part.isdigit():
                return int(chr_part)
            elif chr_part == "X":
                return int(23)  # Put ChrX at the end
            elif chr_part == "Y":
                return int(24)  # Put ChrY after ChrX
            return key

        class GenomeIterator:
            def __init__(self, genome: Genome):
                self.genome = genome
                self.KT_slots = genome.full_KT
                self.KT_slot_keys = sorted(genome.full_KT.keys(), key=custom_sort_chr)
                self.current_slot_index = 0
                self.current_chromosome_index = 0

            def __next__(self):
                if self.current_slot_index < len(self.KT_slots):
                    current_slot = self.KT_slots[self.KT_slot_keys[self.current_slot_index]]
                    if self.current_chromosome_index < len(current_slot):
                        chromosome = current_slot[self.current_chromosome_index]
                        self.current_chromosome_index += 1
                        return chromosome
                    else:
                        self.current_slot_index += 1
                        self.current_chromosome_index = 0
                        return next(self)
                else:
                    raise StopIteration
        return GenomeIterator(self)

    def get_chromosome_list(self):
        chr_list = []
        for chromosome in self:
            chr_list.append(chromosome)
        return chr_list

    def append_history(self, event_type: str, segments: [Segment], chr_from: Chromosome, chr_to: Chromosome):
        """
        add the newest event to history log
        :param event_type: full name of event
        :param segments: segments that were selected for the event (e.g. if +8 inverted to -8, +8 is recorded)
        :param chr_to: object reference
        :param chr_from: object reference
        :return: None
        """
        new_history = tuple([event_type, Arm(segments), chr_from, chr_to])
        self.history.append(new_history)

    def mark_history(self, block_name):
        last_event_in_block = len(self.history) - 1
        self.history_block_markings[last_event_in_block] = block_name

    def history_tostring(self):
        segment_dict = self.segment_indexing()
        return_str = self.initialization_string + '\n'

        sorted_block_marking_keys = sorted(self.history_block_markings)
        start_index = 0
        block_counter = 1
        for history_block_itr in sorted_block_marking_keys:
            return_str += 'block {}: {}\n'.format(str(block_counter), self.history_block_markings[history_block_itr])
            for history_itr_index in range(start_index, history_block_itr + 1):
                segment_str = []
                history_itr = self.history[history_itr_index]
                for segment_itr in history_itr[1].segments:
                    if segment_itr.direction():
                        segment_str.append(segment_dict[segment_itr] + '+')
                    else:
                        new_segment_itr = segment_itr.duplicate()
                        new_segment_itr.invert()
                        segment_str.append(segment_dict[new_segment_itr] + '+')
                return_str += '\t{} on segments [{}], from {} to {}\n'.format(history_itr[0], ','.join(segment_str),
                                                                              history_itr[2].name, history_itr[3].name)
            start_index = history_block_itr + 1
            block_counter += 1
        return return_str

    def segment_indexing(self):
        segment_dict = {}
        current_index = 1
        for segment_itr in self.motherboard.segments:
            segment_dict[segment_itr] = str(current_index)
            current_index += 1
        for centromere_itr in self.centromere_segments:
            centromere_name = centromere_itr.chr_name.replace('Chr', 'CEN')
            segment_dict[centromere_itr] = centromere_name
        return segment_dict

    def motherboard_tostring(self):
        segment_dict = self.segment_indexing()
        sorted_segments = sorted(segment_dict)
        return_str = 'index\torigin\tstart\tend\n'
        for segment_itr in sorted_segments:
            return_str += '{}\t{}\t{}\t{}\n'.format(segment_dict[segment_itr], segment_itr.chr_name,
                                                    segment_itr.start, segment_itr.end)
        return return_str

    def KT_tostring(self):
        segment_dict = self.segment_indexing()
        return_str = 'chromosome\tKT\ttelo1_len\ttelo2_len\n'

        for chr_itr in self:
            if chr_itr.deleted:
                return_str += '{}\tdeleted\t0\t0\n'.format(chr_itr.name)
                continue

            tostring_segment_list = []
            for segment_itr in chr_itr:
                if segment_itr.direction():
                    tostring_segment_list.append(segment_dict[segment_itr] + '+')
                else:
                    new_segment_itr = segment_itr.duplicate()
                    new_segment_itr.invert()
                    tostring_segment_list.append(segment_dict[new_segment_itr] + '-')

            return_str += '{}\t{}\t{}\t{}\n'.format(chr_itr.name, ','.join(tostring_segment_list),
                                                    str(chr_itr.t1_len), str(chr_itr.t2_len))
        return return_str

    def generate_breakpoint(self, event_arm: Arm, breakpoint_index: int):
        """
        split segment such that the breakpoint_index is garenteed to be the end index of a Segment
        :param event_arm: Arm which the event happens on, and the breakpoint_index point at
        :param breakpoint_index: the position of break on the current Arm
            (left_event_index - 1) OR (right_event_index)
        :return: None
        """
        if breakpoint_index == -1:
            # this happens when the break point is at the very beginning of the event_arm, no breaking required
            return

        segment_to_break = Segment('temp', -1, -1)
        left_delete_len = -1
        right_delete_len = -1

        current_bp_index = -1  # corrects 0-index off-shift

        # locate the Segment to create breakpoint
        for segment in event_arm.segments:
            current_bp_index += len(segment)
            if current_bp_index == breakpoint_index:
                # breakpoint exists
                return
            elif current_bp_index > breakpoint_index:
                # document the breakpoint location on the current segment
                segment_to_break = segment.duplicate()
                previous_bp_index = current_bp_index - len(segment)
                left_delete_len = breakpoint_index - previous_bp_index
                right_delete_len = current_bp_index - breakpoint_index
                break
            else:
                # breakpoint location not yet met
                continue

        # create breakpoint on all identical Segments in the genome
        def break_segment(current_arm: Arm):
            left_segment = segment_to_break.duplicate()
            right_segment = segment_to_break.duplicate()
            left_segment.right_delete(right_delete_len)
            right_segment.left_delete(left_delete_len)

            same_direction_match = \
                [index for index, value in enumerate(current_arm.segments) if value == segment_to_break]
            for segment_index_itr in reversed(same_direction_match):
                current_arm.segments.pop(segment_index_itr)
                current_arm.segments.insert(segment_index_itr, left_segment.duplicate())
                current_arm.segments.insert(segment_index_itr + 1, right_segment.duplicate())

            all_match = \
                [index for index, value in enumerate(current_arm.segments)
                 if segment_to_break.same_segment_ignore_dir(value)]
            same_direction_match = \
                [index for index, value in enumerate(current_arm.segments) if value == segment_to_break]
            reversed_direction_match = [element for element in all_match if element not in same_direction_match]
            for segment_index_itr in reversed(reversed_direction_match):
                current_arm.segments.pop(segment_index_itr)
                new_right_segment = right_segment.duplicate()
                new_left_segment = left_segment.duplicate()
                new_right_segment.invert()
                new_left_segment.invert()
                current_arm.segments.insert(segment_index_itr, new_right_segment)
                current_arm.segments.insert(segment_index_itr + 1, new_left_segment)

        break_segment(self.motherboard)
        for slot in self.full_KT:
            for chromosome in self.full_KT[slot]:
                break_segment(chromosome.p_arm)
                break_segment(chromosome.q_arm)
        for history_itr in self.history:
            event_arm_itr = history_itr[1]
            break_segment(event_arm_itr)

    def locate_segments_for_event(self, event_arm: Arm, left_event_index: int, right_event_index: int) \
            -> ([Segment], [int]):
        """
        create breakpoint and select the Segment between the breakpoints
        :param event_arm: chromosome arm that the event will happen in
        :param left_event_index: beginning of deletion, this index will be deleted
        :param right_event_index: end of deletion, this index will be deleted; if -1, then only left_event_index is used
        :return: a list of Segment for processing the event
        """
        segments_selected = []
        segments_selected_indices = []

        if right_event_index >= 0:
            # for an interval
            self.generate_breakpoint(event_arm, left_event_index - 1)
            self.generate_breakpoint(event_arm, right_event_index)

            current_segment_index = 0
            current_bp_index = -1  # corrects 0-index off-shift
            for segment in event_arm.segments:
                current_bp_index += len(segment)
                if left_event_index <= current_bp_index <= right_event_index:
                    segments_selected.append(segment)
                    segments_selected_indices.append(current_segment_index)
                elif current_bp_index > right_event_index:
                    break
                current_segment_index += 1
        elif right_event_index == -1:
            # for a single index, likely as an insertion location
            self.generate_breakpoint(event_arm, left_event_index - 1)

            current_segment_index = 0
            current_bp_index = -1  # corrects 0-index off-shift
            for segment in event_arm.segments:
                current_bp_index += len(segment)
                if left_event_index <= current_bp_index:
                    segments_selected.append(segment)
                    segments_selected_indices.append(current_segment_index)
                    break
                current_segment_index += 1
        else:
            raise ValueError('locate segment for event index right_event_index error')

        return segments_selected, segments_selected_indices

    def deletion(self, event_chromosome: Chromosome, event_arm: Arm, left_event_index: int, right_event_index: int):
        """
        perform deletion event, inplace, between the two indices (inclusive)
        :param event_chromosome: chromosome that the Arm is located on
        :param event_arm: chromosome arm that the event will happen in
        :param left_event_index: beginning of deletion, this index will be deleted
        :param right_event_index: end of deletion, this index will be deleted
        :return: None
        """
        event_segments, event_segments_indices = \
            self.locate_segments_for_event(event_arm, left_event_index, right_event_index)
        # remove empty segments
        event_arm.delete_segments_by_index(event_segments_indices)
        return event_segments

    def tandem_duplication(self, event_chromosome: Chromosome, event_arm: Arm, left_event_index: int, right_event_index: int):
        """
        duplication even, inplace
        :param event_chromosome: chromosome that the Arm is located on
        :param event_arm: chromosome arm that the event will happen in
        :param left_event_index: beginning of deletion, this index will be deleted
        :param right_event_index: end of deletion, this index will be deleted
        :return: None
        """
        event_segments, event_segment_indices = \
            self.locate_segments_for_event(event_arm, left_event_index, right_event_index)
        # duplicate segments
        event_arm.duplicate_segments_by_index(event_segment_indices)
        return event_segments

    def inversion(self, event_chromosome: Chromosome, event_arm: Arm, left_event_index: int, right_event_index: int):
        """
        inversion even, inplace
        :param event_chromosome: chromosome that the Arm is located on
        :param event_arm: chromosome arm tha t the event will happen in
        :param left_event_index: beginning of deletion, this index will be deleted
        :param right_event_index: end of deletion, this index will be deleted
        :return: None
        """
        event_segments, event_segments_indices = \
            self.locate_segments_for_event(event_arm, left_event_index, right_event_index)
        # invert segments
        event_arm.invert_segments_by_index(event_segments_indices)
        return event_segments

    def right_duplication_inversion(self, event_chromosome: Chromosome, event_arm: Arm,
                                    left_event_index: int, right_event_index: int):
        event_segments, event_segment_indices = self.locate_segments_for_event(event_arm, left_event_index,
                                                                               right_event_index)
        new_segment_start_index = event_segment_indices[-1] + 1
        new_segment_end_index = new_segment_start_index + len(event_segments) - 1
        event_arm.duplicate_segments_by_index(event_segment_indices)
        segments_for_inversion_indices = range(new_segment_start_index, new_segment_end_index + 1)
        event_arm.invert_segments_by_index(segments_for_inversion_indices)
        return event_segments

    def left_duplication_inversion(self, event_chromosome: Chromosome, event_arm: Arm,
                                   left_event_index: int, right_event_index: int):
        event_segments, event_segment_indices = self.locate_segments_for_event(event_arm, left_event_index,
                                                                               right_event_index)
        event_arm.duplicate_segments_by_index(event_segment_indices)
        event_arm.invert_segments_by_index(event_segment_indices)
        return event_segments

    def translocation_reciprocal(self,
                                 event_chromosome1, event_arm1: Arm, arm1_left_index: int, arm1_right_index: int,
                                 event_chromosome2, event_arm2: Arm, arm2_left_index: int, arm2_right_index: int):
        arm1_segments, arm1_segment_indices = \
            self.locate_segments_for_event(event_arm1, arm1_left_index, arm1_right_index)
        arm2_segments, arm2_segment_indices = \
            self.locate_segments_for_event(event_arm2, arm2_left_index, arm2_right_index)
        arm1_start_segment_index = arm1_segment_indices[0]
        arm2_start_segment_index = arm2_segment_indices[0]

        event_arm1.delete_segments_by_index(arm1_segment_indices)
        event_arm2.delete_segments_by_index(arm2_segment_indices)
        event_arm2.segments[arm2_start_segment_index:arm2_start_segment_index] = arm1_segments
        event_arm1.segments[arm1_start_segment_index:arm1_start_segment_index] = arm2_segments
        return [arm1_segments, arm2_segments]

    def translocation_nonreciprocal(self,
                                    event_chromosome1, event_arm1: Arm, arm1_left_index: int, arm1_right_index: int,
                                    event_chromosome2, event_arm2: Arm, arm2_index: int):
        arm1_segments, arm1_segment_indices = \
            self.locate_segments_for_event(event_arm1, arm1_left_index, arm1_right_index)
        arm2_segments, arm2_segment_indices = \
            self.locate_segments_for_event(event_arm2, arm2_index, -1)
        arm2_start_segment_index = arm2_segment_indices[0]

        self.append_history('nonreciprocal translocation', arm1_segments, event_chromosome1, event_chromosome2)

        event_arm1.delete_segments_by_index(arm1_segment_indices)
        event_arm2.segments[arm2_start_segment_index:arm2_start_segment_index] = arm1_segments

    def chromosomal_deletion(self, event_chromosome: Chromosome):
        event_chromosome.deleted = True
        event_segments = \
            event_chromosome.p_arm.segments + event_chromosome.centromere.segments + event_chromosome.q_arm.segments
        return event_segments

    def chromosomal_duplication(self, event_chromosome: Chromosome):
        if event_chromosome.deleted:
            raise ValueError('selected chromosome for duplication is already deleted')
        new_chromosome = event_chromosome.duplicate()
        self.full_KT[new_chromosome.name[:-1]].append(new_chromosome)
        new_chromosome.name = new_chromosome.name[:-1] + chr(len(self.full_KT[new_chromosome.name[:-1]]) + 96)
        event_segments = \
            event_chromosome.p_arm.segments + event_chromosome.centromere.segments + event_chromosome.q_arm.segments
        return event_segments

    def arm_deletion(self, event_chromosome: Chromosome, event_arm: Arm):
        event_segments, event_segments_indices = \
            self.locate_segments_for_event(event_arm, 0, len(event_arm) - 1)
        # remove empty segments
        event_arm.delete_segments_by_index(event_segments_indices)
        return event_segments

    def arm_tandem_duplication(self, event_chromosome: Chromosome, event_arm: Arm):
        event_segments, event_segment_indices = \
            self.locate_segments_for_event(event_arm, 0, len(event_arm) - 1)
        # duplicate segments
        event_arm.duplicate_segments_by_index(event_segment_indices)
        return event_segments

    def output_KT(self, output_file):
        with open(output_file, 'w') as fp_write:
            fp_write.write(self.motherboard_tostring())
            fp_write.write('---\n')
            fp_write.write(self.KT_tostring())
            fp_write.write('---\n')
            fp_write.write(self.history_tostring())

    def output_FASTA(self, genome_path: str, output_file: str):
        """
        :param genome_path: require the header names to match the genome's chromosome names
        :param output_file:
        :return:
        """
        def reverse_complement(dna_sequence):
            complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C',
                               'a': 't', 't': 'a', 'c': 'g', 'g': 'c',
                               'M': 'K', 'm': 'k', 'K': 'M', 'k': 'm',
                               'R': 'Y', 'r': 'y', 'Y': 'R', 'y': 'r',
                               'W': 'W', 'w': 'w',
                               'S': 'S', 's': 's',
                               'B': 'V', 'b': 'v', 'V': 'B', 'v': 'b',
                               'H': 'D', 'h': 'd', 'D': 'H', 'd': 'h',
                               'N': 'N', 'n': 'n',}
            reverse_sequence = dna_sequence[::-1]
            complement_sequence = ''.join(complement_dict[base] for base in reverse_sequence)
            return complement_sequence

        sequence_dict = IO.read_FASTA(genome_path, ['all'])
        output_dict = {}
        for chromosome in self:
            if chromosome.deleted:
                continue

            # telomere 1
            new_sequence = ['N' * chromosome.t1_len]

            # p-arm, centromere, q-arm
            for segment in chromosome:
                if segment.direction():
                    new_sequence.append(sequence_dict[segment.chr_name][segment.start: segment.end + 1])
                else:
                    new_sequence.append(reverse_complement
                                        (sequence_dict[segment.chr_name][segment.end: segment.start + 1]))

            # telomere 2
            new_sequence.append('N' * chromosome.t2_len)

            output_dict[chromosome.name] = ''.join(new_sequence)

        IO.sequence_dict_to_FASTA(output_dict, output_file)
