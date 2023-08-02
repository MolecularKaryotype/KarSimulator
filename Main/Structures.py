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

    def __init__(self, name: str, p_arm: Arm, q_arm: Arm, t1_len: int, t2_len: int, centromere: Arm):
        self.name = name
        self.p_arm = p_arm
        self.q_arm = q_arm
        self.t1_len = t1_len
        self.t2_len = t2_len
        self.centromere = centromere

    def __len__(self):
        return self.p_arm_len() + self.q_arm_len()

    def __str__(self):
        return_str = '{}: t1-{} t2-{}\n\tp-arm: {}\n\tq-arm: {}\n\tCEN: {}' \
            .format(self.name, self.t1_len, self.t2_len, str(self.p_arm), str(self.q_arm), str(self.centromere))
        return return_str

    def p_arm_len(self):
        return len(self.p_arm)

    def q_arm_len(self):
        return len(self.q_arm)

    def duplicate(self):
        return Chromosome(self.name, self.p_arm.duplicate(), self.q_arm.duplicate(),
                          self.t1_len, self.t2_len, self.centromere.duplicate())


class Genome:
    full_KT: {str: [Chromosome]}  # has exactly 24 slots, corresponding to the 24 possible chromosome types
    motherboard: Arm  # using the Arm object to use generate breakpoint method
    centromere_segments = [Segment]
    history: [(str, Arm, Chromosome, Chromosome)]  # event type, event segments, chr from, chr to

    def __init__(self, copy_number: int, chr_of_interest: [str], genome_index_file: str):
        # construct KT slots
        self.full_KT = {}
        self.motherboard = Arm([])
        self.centromere_segments = []
        self.history = []
        for slot in chr_of_interest:
            self.full_KT[slot] = []
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

                    self.full_KT[chr_name].append(Chromosome(chr_name, Arm([p_arm_segment]), Arm([q_arm_segment]),
                                                             t1_len, t2_len, Arm([centromere_segment])))
                    self.motherboard.segments.extend([p_arm_segment.duplicate(),
                                                      q_arm_segment.duplicate()])
                    self.centromere_segments.append(centromere_segment.duplicate())

        # setup copy numbers
        if copy_number <= 0:
            raise ValueError('copy number')
        for slot in self.full_KT:
            first_chromosome = self.full_KT[slot][0]
            for i in range(copy_number - 1):
                self.full_KT[slot].append(first_chromosome.duplicate())
        # rename each chromosome with unique name
        for slot in self.full_KT:
            current_slot = self.full_KT[slot]
            for i in range(copy_number):
                current_slot[i].name = current_slot[i].name + chr(i + 97)

    def __str__(self):
        return_str = ''
        for slot in self.full_KT:
            for chromosome in self.full_KT[slot]:
                return_str += str(chromosome) + '\n'
        return return_str

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

    def history_tostring(self):
        segment_dict = self.segment_indexing()
        return_str = ''
        for history_itr in self.history:
            segment_str = ''
            for segment_itr in history_itr[1].segments:
                if segment_itr.direction():
                    segment_str += segment_dict[segment_itr] + '+'
                else:
                    new_segment_itr = segment_itr.duplicate()
                    new_segment_itr.invert()
                    segment_str += segment_dict[new_segment_itr] + '-'
            return_str += '{} on segments [{}], from {} to {}\n'.format(history_itr[0], segment_str,
                                                                        history_itr[2].name, history_itr[3].name)
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
        for slot in self.full_KT:
            for chr_itr in self.full_KT[slot]:
                tostring_segment_list = []
                for segment_itr in chr_itr.p_arm.segments:
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
        :param right_event_index: end of deletion, this index will be deleted
        :return: a list of Segment for processing the event
        """
        self.generate_breakpoint(event_arm, left_event_index - 1)
        self.generate_breakpoint(event_arm, right_event_index)

        segments_selected = []
        segments_selected_indices = []
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
        # document segments deleted
        self.append_history('deletion', event_segments, event_chromosome, event_chromosome)
        # remove empty segments
        event_arm.delete_segments_by_index(event_segments_indices)

    def duplication(self, event_chromosome: Chromosome, event_arm: Arm, left_event_index: int, right_event_index: int):
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
        # document segments duplicated
        self.append_history('duplication', event_segments, event_chromosome, event_chromosome)
        # duplicate segments
        event_arm.duplicate_segments_by_index(event_segment_indices)

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
        # document segments inverted
        self.append_history('inversion', event_segments, event_chromosome, event_chromosome)
        # invert segments
        event_arm.invert_segments_by_index(event_segments_indices)
