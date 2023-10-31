from sequence_dict_to_FASTA import *
from read_FASTA import read_FASTA
import copy


class Segment:
    chr_name: str
    start: int
    end: int
    segment_type: str
    kt_index: str
    ordinal: int

    def __init__(self, chr_name: str, start: int, end: int, segment_type=None, kt_index=None):
        self.chr_name = chr_name
        self.start = start
        self.end = end
        self.segment_type = segment_type
        self.kt_index = kt_index
        self.ordinal = -1

    def __len__(self):
        return abs(self.end - self.start) + 1

    def __lt__(self, other):
        type_order = ["telomere1", "telomere2", "centromere", "hardmask", "superdup"]

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
            if max(self.start, self.end) != max(other.start, other.end):
                return max(self.start, self.end) < max(other.start, other.end)
            else:
                self_type_index = type_order.index(self.segment_type)
                other_type_index = type_order.index(other.segment_type)
                return self_type_index < other_type_index

    def __eq__(self, other):
        if isinstance(other, Segment):
            return (self.chr_name, self.start, self.end) == (other.chr_name, other.start, other.end)
        return False
        # return (self.chr_name, self.start, self.end) == (other.chr_name, other.start, other.end)

    def __hash__(self):
        return hash((self.chr_name, self.start, self.end))

    def __str__(self):
        additional_info = ""
        if self.segment_type is not None:
            additional_info += ", " + self.segment_type
        if self.kt_index is not None:
            additional_info += ", " + self.kt_index
        if self.ordinal != -1:
            additional_info += ", " + str(self.ordinal)
        return_str = "({}, {}, {}{})".format(self.chr_name, self.start, self.end, additional_info)
        return return_str

    def annotated_number(self):
        start_str = "{:,}".format(self.start)
        end_str = "{:,}".format(self.end)
        return "({}-{}-{}-{})".format(self.chr_name, start_str, end_str, self.segment_type)

    def alignment_output(self):
        start_str = "{:,}".format(self.start)
        end_str = "{:,}".format(self.end)
        if self.segment_type is None and self.kt_index is None:
            return "({} {} {})".format(self.chr_name, start_str, end_str)
        elif self.segment_type is None:
            return "({} {} {} {})".format(self.chr_name, start_str, end_str, self.kt_index)
        else:
            return "({} {} {} {} {})".format(self.chr_name, start_str, end_str, self.segment_type, self.kt_index)

    def concise_str(self):
        if self.segment_type is None:
            s1 = ""
        else:
            s1 = self.segment_type
        if self.kt_index is None:
            s2 = ""
        else:
            s2 = self.kt_index
        return "({},{})".format(s1, s2)

    def thousand_delimited(self):
        return "({}-{}-{})".format(self.chr_name, "{:,}".format(self.start), "{:,}".format(self.end))

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
        return Segment(self.chr_name, self.start, self.end, self.segment_type, self.kt_index)

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

    def invert(self, inplace=True):
        if inplace:
            temp_start = self.start
            self.start = self.end
            self.end = temp_start
        else:
            new_segment = self.duplicate()
            new_segment.start = self.end
            new_segment.end = self.start
            return new_segment

    def segment_intersection(self, other_segment):
        duplicate_self = self.duplicate()
        duplicate_other = other_segment.duplicate()
        if not duplicate_self.direction():
            duplicate_self.invert()
        if not duplicate_other.direction():
            duplicate_other.invert()

        if duplicate_self.chr_name != duplicate_other.chr_name:
            return False
        if duplicate_self.start <= duplicate_other.end and duplicate_self.end >= duplicate_other.start:
            return True
        else:
            return False

    def bp_in_interior(self, bp_chromosome, bp_index, bp_type):
        """
        For KarComparator
        Used to see if an index is within the segment, potentially requiring a breakpoint
        :param bp_chromosome:
        :param bp_index:
        :param bp_type: start or end
        :return:
        """
        if self.chr_name != bp_chromosome:
            return False
        if self.direction():
            if bp_type == "start" and bp_index == self.start:
                return False
            elif bp_type == "end" and bp_index == self.end:
                return False

            if bp_index in range(self.start, self.end + 1):
                return True
            else:
                return False
        else:
            if bp_type == "end" and bp_index == self.start:
                return False
            elif bp_type == "start" and bp_index == self.end:
                return False

            if bp_index in range(self.end, self.start + 1):
                return True
            else:
                return False


class Arm:
    segments: [Segment]
    deleted: bool
    arm_type: str

    def __init__(self, segments: [Segment], arm_type: str):
        self.segments = segments
        self.deleted = False
        self.arm_type = arm_type

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

    def get_segment_indices(self, input_data):
        """
        :param input_data: Segment or [Segment]
        :return: int or [int]
        """
        # the same segment ID (segment object) will occur only once per Genome, no exception
        if isinstance(input_data, Segment):
            for arm_segment_ind in range(0, len(self.segments)):
                if self.segments[arm_segment_ind] is input_data:
                    return arm_segment_ind
        elif isinstance(input_data, list) and all(isinstance(item, Segment) for item in input_data):
            return_indices = []
            for segment_to_find in input_data:
                for arm_segment_ind in range(0, len(self.segments)):
                    if self.segments[arm_segment_ind] is segment_to_find:
                        return_indices.append(arm_segment_ind)
                        break
            return return_indices
        else:
            raise TypeError('get_segment_indices wrong type')

    def duplicate(self):
        new_segments = []
        for segment in self.segments:
            new_segments.append(segment.duplicate())
        return Arm(new_segments, self.arm_type)

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

    def arm_intersection(self, other_arm):
        for segment1 in self.segments:
            for segment2 in other_arm.segments:
                if segment1.segment_intersection(segment2):
                    return True
        return False

    def report_arm_intersection(self, other_arm):
        intersecting_segments = []
        for segment1 in self.segments:
            for segment2 in other_arm.segments:
                if segment1.segment_intersection(segment2):
                    intersecting_segments.append(segment2)
        return_str = ''
        for segment in intersecting_segments:
            return_str += segment.annotated_number()
        return return_str

    def gather_boundary_points(self):
        """
        For KarComparator
        :return: a list of all the boundary points for each Segment (two for each)
        """
        return_list = []
        for segment in self.segments:
            if segment.direction():
                return_list.append(tuple([segment.chr_name, segment.start, 'start']))
                return_list.append(tuple([segment.chr_name, segment.end, 'end']))
            else:
                return_list.append(tuple([segment.chr_name, segment.start, 'end']))
                return_list.append(tuple([segment.chr_name, segment.end, 'start']))
        return return_list

    def introduce_breakpoint(self, bp_chromosome, bp_index, bp_type):
        """
        For KarComparator
        Search through the arm and generate the breakpoint, if within an interior of a Segment
        :param bp_chromosome:
        :param bp_index:
        :param bp_type: start or end bp
        :return:
        """
        current_segment_index = 0
        while current_segment_index < len(self.segments):
            current_segment = self.segments[current_segment_index]
            if current_segment.bp_in_interior(bp_chromosome, bp_index, bp_type):
                insertion_index = self.get_segment_index(current_segment)
                if current_segment.direction():
                    if bp_type == "start":
                        left_segment = Segment(current_segment.chr_name, current_segment.start, bp_index - 1,
                                               current_segment.segment_type, current_segment.kt_index)
                        right_segment = Segment(current_segment.chr_name, bp_index, current_segment.end,
                                                current_segment.segment_type, current_segment.kt_index)
                    elif bp_type == "end":
                        left_segment = Segment(current_segment.chr_name, current_segment.start, bp_index,
                                               current_segment.segment_type, current_segment.kt_index)
                        right_segment = Segment(current_segment.chr_name, bp_index + 1, current_segment.end,
                                                current_segment.segment_type, current_segment.kt_index)
                    else:
                        raise ValueError('bp_type must be start OR end')
                else:
                    if bp_type == "start":
                        left_segment = Segment(current_segment.chr_name, current_segment.start, bp_index,
                                               current_segment.segment_type, current_segment.kt_index)
                        right_segment = Segment(current_segment.chr_name, bp_index - 1, current_segment.end,
                                                current_segment.segment_type, current_segment.kt_index)
                    elif bp_type == "end":
                        left_segment = Segment(current_segment.chr_name, current_segment.start, bp_index + 1,
                                               current_segment.segment_type, current_segment.kt_index)
                        right_segment = Segment(current_segment.chr_name, bp_index, current_segment.end,
                                                current_segment.segment_type, current_segment.kt_index)
                    else:
                        raise ValueError('bp_type must be start OR end')

                self.segments.pop(insertion_index)
                self.segments.insert(insertion_index, right_segment)
                self.segments.insert(insertion_index, left_segment)
                current_segment_index += 1  # since we added one more segment in-place

            current_segment_index += 1

    def get_segment_index(self, input_segment):
        """
        For KarComparator
        find the index of segment in list, matching the object using __is__ (not __eq__)
        :param input_segment: Only use segment that is in the Arm
        :return:
        """
        for segment_index in range(0, len(self.segments)):
            current_segment = self.segments[segment_index]
            if current_segment is input_segment:
                return segment_index
        raise RuntimeError('segment not found in Arm')

    def merge_breakpoints(self):
        """
        For Masking File Generation
        :return:
        """
        current_segment_index = 0
        while current_segment_index < len(self.segments) - 1:
            current_segment = self.segments[current_segment_index]
            next_segment = self.segments[current_segment_index + 1]
            if current_segment.chr_name == next_segment.chr_name:
                if current_segment.segment_type == next_segment.segment_type:
                    if current_segment.direction() and current_segment.end + 1 == next_segment.start:
                        new_segment = Segment(current_segment.chr_name, current_segment.start,
                                              next_segment.end, current_segment.segment_type)
                        self.segments.pop(current_segment_index)
                        self.segments.pop(current_segment_index)
                        self.segments.insert(current_segment_index, new_segment)
                        continue
                    elif (not current_segment.direction()) and current_segment.end - 1 == next_segment.start:
                        new_segment = Segment(current_segment.chr_name, current_segment.start,
                                              next_segment.end, current_segment.segment_type)
                        self.segments.pop(current_segment_index)
                        self.segments.pop(current_segment_index)
                        self.segments.insert(current_segment_index, new_segment)
                        continue

            current_segment_index += 1


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
        if self.deleted:
            return 0
        else:
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
                    if current_arm.deleted:
                        self.current_arm_index += 1
                        return next(self)
                    elif self.current_segment_index < len(current_arm.segments):
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
        # TODO: implement
        pass

    def p_arm_len(self):
        if self.p_arm.deleted:
            return 0
        else:
            return len(self.p_arm)

    def q_arm_len(self):
        if self.q_arm.deleted:
            return 0
        else:
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
    pass  # TODO: implement


class Genome:
    full_KT: {str: [Chromosome]}  # has as many slots as there are chromosome type, i.e. 24 for a male, 23 for a female
    motherboard: Arm  # using the Arm object to use generate breakpoint method
    centromere_segments = [Segment]
    history_block_markings = {}  # history enumerate index: block name
    history: [(str, Arm, Chromosome, Chromosome)]  # event type, event segments, chr from, chr to
    ordinal_history: [[(str, Segment, Segment, Segment)]]  # ins/del/inv, event_segment, left_bound, right_bound
    initialization_string: str  # contains information with the initialization of the genome

    def __init__(self, full_KT, motherboard_segments, centromere_segments, initialization_string,
                 history=None, history_markers=None, ordinal_history=None):
        self.full_KT = full_KT
        self.motherboard = Arm(motherboard_segments, 'motherboard')
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
        if ordinal_history is not None:
            self.ordinal_history = ordinal_history
        else:
            self.ordinal_history = []

    def duplicate(self):
        new_full_KT = {}
        for key in self.full_KT:
            new_chr_list = []
            for chr_itr in self.full_KT[key]:
                new_chr_list.append(chr_itr.duplicate())
            new_full_KT[key] = new_chr_list

        new_centromere_segments = []
        for chr_itr in self.centromere_segments:
            new_centromere_segments.append(chr_itr.duplicate())

        new_history = []
        for history_itr in self.history:
            new_item = [copy.deepcopy(history_itr[0])]

        return Genome(new_full_KT,
                      self.motherboard.duplicate(),
                      new_centromere_segments,
                      copy.deepcopy(self.initialization_string),
                      )

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
        new_history = tuple([event_type, Arm(segments, 'history'), chr_from, chr_to])
        self.history.append(new_history)

    def mark_history(self, block_name):
        last_event_in_block = len(self.history) - 1
        self.history_block_markings[last_event_in_block] = block_name

    def pop_last_history_marking(self):
        last_event_in_block = len(self.history) - 1
        self.history_block_markings.pop(last_event_in_block)

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

    def need_breakpoint(self, event_arm: Arm, breakpoint_index: int):
        """
        split segment such that the breakpoint_index is guarenteed to be the end index of a Segment
        :param event_arm: Arm which the event happens on, and the breakpoint_index point at
        :param breakpoint_index: the position of break on the current Arm
            (left_event_index - 1) OR (right_event_index)
        :return: None
        """
        if breakpoint_index == -1:
            # this happens when the break point is at the very beginning of the event_arm, no breaking required
            return 0

        segment_to_break = Segment('temp', -1, -1)
        left_delete_len = -1
        right_delete_len = -1

        current_bp_index = -1  # corrects 0-index off-shift

        # locate the Segment to create breakpoint
        for segment in event_arm.segments:
            current_bp_index += len(segment)
            if current_bp_index == breakpoint_index:
                # breakpoint exists
                return 0
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
        return 1

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
        :param left_event_index: beginning of event, this index will be included
        :param right_event_index: end of event, this index will be included; if -1, then only left_event_index is used
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

    def deletion(self, event_arm: Arm, left_event_index: int, right_event_index: int):
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

    def tandem_duplication(self, event_arm: Arm, left_event_index: int, right_event_index: int):
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

    def inversion(self, event_arm: Arm, left_event_index: int, right_event_index: int):
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

    def right_duplication_inversion(self, event_arm: Arm, left_event_index: int, right_event_index: int):
        event_segments, event_segment_indices = self.locate_segments_for_event(event_arm, left_event_index,
                                                                               right_event_index)
        new_segment_start_index = event_segment_indices[-1] + 1
        new_segment_end_index = new_segment_start_index + len(event_segments) - 1
        event_arm.duplicate_segments_by_index(event_segment_indices)
        segments_for_inversion_indices = range(new_segment_start_index, new_segment_end_index + 1)
        event_arm.invert_segments_by_index(segments_for_inversion_indices)
        return event_segments

    def left_duplication_inversion(self, event_arm: Arm,
                                   left_event_index: int, right_event_index: int):
        event_segments, event_segment_indices = self.locate_segments_for_event(event_arm, left_event_index,
                                                                               right_event_index)
        event_arm.duplicate_segments_by_index(event_segment_indices)
        event_arm.invert_segments_by_index(event_segment_indices)
        return event_segments

    def translocation_reciprocal_balanced(self,
                                          event_arm1: Arm, arm1_left_index: int, arm1_right_index: int,
                                          event_arm2: Arm, arm2_left_index: int, arm2_right_index: int):
        arm1_segments, _ = \
            self.locate_segments_for_event(event_arm1, arm1_left_index, arm1_right_index)
        arm2_segments, _ = \
            self.locate_segments_for_event(event_arm2, arm2_left_index, arm2_right_index)

        # notate the segment right before our segments, None if our segment is the first
        # this fix the problem when two segments are on the same arm, creating index-shifts
        # this is reliable because
        # 1) segment object ID is ALWAYS unique on a genome
        # 2) the range of segments are always continuous
        arm1_start_segment_index = event_arm1.get_segment_indices(arm1_segments[0])
        arm2_start_segment_index = event_arm2.get_segment_indices(arm2_segments[0])
        arm1_prior_segment = 'None'
        arm2_prior_segment = 'None'
        if arm1_start_segment_index != 0:
            arm1_prior_segment = event_arm1.segments[arm1_start_segment_index - 1]
        if arm2_start_segment_index != 0:
            arm2_prior_segment = event_arm2.segments[arm2_start_segment_index - 1]

        arm1_segment_indices = event_arm1.get_segment_indices(arm1_segments)
        event_arm1.delete_segments_by_index(arm1_segment_indices)
        arm2_segment_indices = event_arm2.get_segment_indices(arm2_segments)
        event_arm2.delete_segments_by_index(arm2_segment_indices)

        if arm1_prior_segment == 'None':
            arm1_add_index = 0
        else:
            arm1_add_index = event_arm1.get_segment_indices(arm1_prior_segment) + 1
        if arm2_prior_segment == 'None':
            arm2_add_index = 0
        else:
            arm2_add_index = event_arm2.get_segment_indices(arm2_prior_segment) + 1

        event_arm1.segments[arm1_add_index:arm1_add_index] = arm2_segments
        event_arm2.segments[arm2_add_index:arm2_add_index] = arm1_segments

        return [arm1_segments, arm2_segments]

    def translocation_reciprocal_unbalanced(self,
                                            event_arm1: Arm, arm1_left_index: int, arm1_right_index: int,
                                            event_arm2: Arm, arm2_left_index: int, arm2_right_index: int):
        arm1_segments, _ = \
            self.locate_segments_for_event(event_arm1, arm1_left_index, arm1_right_index)
        arm2_segments, _ = \
            self.locate_segments_for_event(event_arm2, arm2_left_index, arm2_right_index)

        arm2_prior_segment = 'None'
        arm2_start_segment_index = event_arm2.get_segment_indices(arm2_segments[0])
        if arm2_start_segment_index != 0:
            arm2_prior_segment = event_arm2.segments[arm2_start_segment_index - 1]

        arm1_segment_indices = event_arm1.get_segment_indices(arm1_segments)
        event_arm1.delete_segments_by_index(arm1_segment_indices)
        arm2_segment_indices = event_arm2.get_segment_indices(arm2_segments)
        event_arm2.delete_segments_by_index(arm2_segment_indices)

        if arm2_prior_segment == 'None':
            arm2_add_index = 0
        else:
            arm2_add_index = event_arm2.get_segment_indices(arm2_prior_segment) + 1

        event_arm2.segments[arm2_add_index:arm2_add_index] = arm1_segments

        return [arm1_segments, arm2_segments]

    def translocation_nonreciprocal(self,
                                    event_arm1: Arm, arm1_left_index: int, arm1_right_index: int,
                                    event_arm2: Arm, arm2_index: int):
        arm1_segments, _ = \
            self.locate_segments_for_event(event_arm1, arm1_left_index, arm1_right_index)
        arm2_latter_segment, _ = \
            self.locate_segments_for_event(event_arm2, arm2_index, -1)

        arm1_segment_indices = event_arm1.get_segment_indices(arm1_segments)
        event_arm1.delete_segments_by_index(arm1_segment_indices)
        arm2_add_index = event_arm2.get_segment_indices(arm2_latter_segment[0])
        event_arm2.segments[arm2_add_index:arm2_add_index] = arm1_segments

        return arm1_segments

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
        event_segments, _ = \
            self.locate_segments_for_event(event_arm, 0, len(event_arm) - 1)
        event_arm.deleted = True
        if len(event_chromosome) == 0:
            event_chromosome.deleted = True
        return event_segments

    def arm_tandem_duplication(self, event_chromosome: Chromosome, event_arm: Arm):
        new_chromosome = event_chromosome.duplicate()
        self.full_KT[new_chromosome.name[:-1]].append(new_chromosome)
        new_chromosome.name = new_chromosome.name[:-1] + chr(len(self.full_KT[new_chromosome.name[:-1]]) + 96)
        event_segments, event_segment_indices = \
            self.locate_segments_for_event(event_arm, 0, len(event_arm) - 1)
        # duplicate arm onto a new chromosome
        if event_arm.arm_type == 'p':
            new_chromosome.q_arm.deleted = True
        else:
            new_chromosome.p_arm.deleted = True
        return event_segments, new_chromosome

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
                               'N': 'N', 'n': 'n', }
            reverse_sequence = dna_sequence[::-1]
            complement_sequence = ''.join(complement_dict[base] for base in reverse_sequence)
            return complement_sequence

        sequence_dict = read_FASTA(genome_path, ['all'])
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

        sequence_dict_to_FASTA(output_dict, output_file)


class Path:
    linear_path: Arm
    path_chr: str
    path_name: str

    def __init__(self, linear_path: Arm, path_name=None, path_chr=None):
        self.linear_path = linear_path
        self.path_chr = path_chr
        self.path_name = path_name

    def __str__(self):
        return str("chr_bin: {}, path_name: {}, segments: {}".format(self.path_chr, self.path_name, self.linear_path))

    def concise_str(self):
        segment_str = ""
        for segment in self.linear_path.segments:
            segment_str += segment.concise_str()
        return str("path_name: {}, segments: {}".format(self.path_name, segment_str))

    def reverse(self):
        new_segments = []
        for segment in reversed(self.linear_path.segments):
            new_segment = segment.duplicate()
            new_segment.invert()
            new_segments.append(new_segment)
        self.linear_path.segments = new_segments

    def get_path_notes(self):
        segment_origin_str = ""
        for segment in self.linear_path.segments:
            segment_origin_str += segment.chr_name + " "
        return str("chr_bin: {}, path_name: {}, segment_chr: {}".format(self.path_chr, self.path_name,
                                                                        segment_origin_str))

    def generate_mutual_breakpoints(self, other_path=None, mutual=True):
        """
        make sure all segments within the 1/2 paths have mutual breakpoints
        :param other_path: if None, then breaking within itself
        :param mutual: whether to generate breakpoints on the other_path
        :return:
        """
        path1_breakpoints = self.linear_path.gather_boundary_points()

        if other_path is not None:
            path2_breakpoints = other_path.linear_path.gather_boundary_points()
            for breakpoint_itr in path2_breakpoints:
                self.linear_path.introduce_breakpoint(*breakpoint_itr)

            if mutual:
                for breakpoint_itr in path1_breakpoints:
                    other_path.linear_path.introduce_breakpoint(*breakpoint_itr)
        else:
            for breakpoint_itr in path1_breakpoints:
                self.linear_path.introduce_breakpoint(*breakpoint_itr)

    def duplicate(self):
        new_arm = self.linear_path.duplicate()
        return Path(new_arm, self.path_chr, self.path_name)


def segment_intersection_test():
    segment1 = Segment("Chr1", 10, 30)
    segment2 = Segment("Chr1", 25, 40)
    segment3 = Segment("Chr2", 5, 15)
    segment4 = Segment("Chr1", 5, 9)
    print(segment1.segment_intersection(segment2))  # Should return True
    print(segment1.segment_intersection(segment3))  # Should return False
    print(segment1.segment_intersection(segment4))  # Should return False


def arm_intersection_test():
    pass


def break_masking_file():
    with open("../Metadata/merged_masking_annotated.bed") as fp_read:
        fp_read.readline()
        segments = []
        for line in fp_read:
            line = line.replace('\n', '').split('\t')
            if line[0] == "23":
                chr_name = "X"
            elif line[0] == "24":
                chr_name = "Y"
            else:
                chr_name = line[0]
            chr_name = "Chr" + chr_name
            segments.append(Segment(chr_name, int(line[1]), int(line[2]), line[3]))

    masking_arm = Arm(segments, "masking_arm")
    masking_path = Path(masking_arm)

    masking_path.generate_mutual_breakpoints()
    masking_path.linear_path.segments = sorted(masking_path.linear_path.segments)

    # remove duplicates
    unique_segments = []
    prev_segment = None

    for segment in masking_path.linear_path.segments:
        if segment != prev_segment:
            unique_segments.append(segment)  # save the first occurrence, assume list is sorted
            prev_segment = segment

    masking_path.linear_path.segments = unique_segments
    masking_path.linear_path.merge_breakpoints()
    print(len(masking_path.linear_path))

    centromere_telomere_length = 0
    for segment in masking_path.linear_path.segments:
        if segment.segment_type in ["centromere", "telomere1", "telomere2"]:
            centromere_telomere_length += len(segment)
    print(centromere_telomere_length)

    hardmask_length = 0
    for segment in masking_path.linear_path.segments:
        if segment.segment_type in ["hardmask"]:
            hardmask_length += len(segment)
    print(hardmask_length)

    superdup_length = 0
    for segment in masking_path.linear_path.segments:
        if segment.segment_type in ["superdup"]:
            superdup_length += len(segment)
    print(superdup_length)

    with open("../Metadata/merged_masking_unique.bed", "w") as fp_write:
        fp_write.write("Chr\tStartPos\tEndPos\tType\n")
        for segment in masking_path.linear_path.segments:
            fp_write.write("{}\t{}\t{}\t{}\n".format(segment.chr_name, segment.start, segment.end,
                                                     segment.segment_type))


if __name__ == "__main__":
    segment_intersection_test()
