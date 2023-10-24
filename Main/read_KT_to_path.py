from Structures import *
from read_masking_regions import read_masking_regions
from Start_Genome import *


def read_KT_to_path(KT_file, masking_file):
    genome = generate_genome_from_KT(KT_file)
    history_dict = extract_history_by_chr_destination(genome)
    index_dict = genome.segment_indexing()
    t2_segments = get_t2_segments(masking_file)
    path_list = []

    history_counter = 1
    for chromosome_itr in genome:
        t1_segment = Segment(chromosome_itr.name[:-1], 0, chromosome_itr.t1_len - 1, "telomere1")
        segment_list = [t1_segment]

        for segment_itr in chromosome_itr.p_arm.segments:
            segment_itr.kt_index = get_kt_index(index_dict, segment_itr)
            if segment_itr.direction():
                segment_itr.kt_index += "+"
            else:
                segment_itr.kt_index += "-"
            segment_list.append(segment_itr)

        cen_segment = chromosome_itr.centromere.segments[0].duplicate()
        cen_segment.segment_type = "centromere"
        cen_segment.kt_index = get_kt_index(index_dict, cen_segment)
        segment_list.append(cen_segment)

        for segment_itr in chromosome_itr.q_arm.segments:
            segment_itr.kt_index = get_kt_index(index_dict, segment_itr)
            if segment_itr.direction():
                segment_itr.kt_index += "+"
            else:
                segment_itr.kt_index += "-"
            segment_list.append(segment_itr)

        t2_segment = None
        for find_t2_segment_itr in t2_segments.segments:
            if find_t2_segment_itr.chr_name == chromosome_itr.name[:-1]:
                t2_segment = find_t2_segment_itr
                break
        segment_list.append(t2_segment)

        # # if no event happened on this chr
        # if chromosome_itr.name not in history_dict:
        #     path_list.append(Path(Arm(segment_list, "KT_path"), chromosome_itr.name, chromosome_itr.name[:-1]))
        #     continue

        # label SV segments from history
        # get inv/ins first
        for ordinal_history_entry_index in range(len(genome.ordinal_history)):
            ordinal_history_entry = genome.ordinal_history[ordinal_history_entry_index]
            for ordinal_sub_entry in ordinal_history_entry:
                if ordinal_sub_entry[0] in ['inv', 'ins']:
                    event_segment = ordinal_sub_entry[1]
                    event_segment_in_path = get_segment_from_ordinal(event_segment, event_segment.ordinal, segment_list)
                    event_segment_in_path.segment_type = ordinal_sub_entry[0] + ": SV" + str(ordinal_history_entry_index)

        # intersperse del in-between the left and right boundary

        path_list.append(Path(Arm(segment_list, "KT_path"), chromosome_itr.name, chromosome_itr.name[:-1]))
        print('x')

    return path_list


def get_kt_index(input_index_dict, input_segment):
    for key in input_index_dict:
        if input_segment.same_segment_ignore_dir(key):
            return input_index_dict[key]


def extract_history_by_chr_destination(genome):
    history_dict = {}
    for history_entry in genome.history:
        history_type = history_entry[0]
        history_arm = history_entry[1]
        history_destination_chr = history_entry[3].name
        if history_destination_chr not in history_dict:
            history_dict[history_destination_chr] = []
        history_dict[history_destination_chr].append(tuple([history_type, history_arm]))
    return history_dict


def get_segment_in_SV(history_entry, segment_list):
    event_name = history_entry[0]
    if event_name == 'deletion':
        return
    # for
    # for segment in segment_list:


def get_t2_segments(masking_file):
    masking_arm = read_masking_regions(masking_file)
    t2_segments = []
    for segment in masking_arm.segments:
        if segment.segment_type == "telomere2":
            t2_segments.append(segment)
    return Arm(t2_segments, "t2_segments")


def segments_are_continuous(segment1: Segment, segment2: Segment):
    if segment1.direction():
        if segment1.end + 1 == segment2.start:
            return True
    else:
        if segment1.end - 1 == segment2.start:
            return True
    return False


def get_segment_from_ordinal(input_current_segment: Segment, input_segment_ordinal: int, input_segment_list: []):
    output_segment = None
    finder_ptr = 0
    p_arm_exhausted = False
    for occurrence in range(0, input_segment_ordinal):
        # FIXME: runtime error require investigation
        segment_not_matched = True
        while segment_not_matched:
            if input_segment_list[finder_ptr].same_segment_ignore_dir(input_current_segment):
                output_segment = input_segment_list[finder_ptr]
                segment_not_matched = False
            finder_ptr += 1
    return output_segment


def test():
    return_list = read_KT_to_path("/media/zhaoyang-new/workspace/KarSim/KarSimulator/test_folder/23Xe10_r1.kt.txt",
                                  "../Metadata/merged_masking_unique.bed")
    for path in return_list:
        print(path.concise_str())


if __name__ == "__main__":
    test()
