from Structures import *
from read_masking_regions import read_masking_regions
from Start_Genome import generate_genome_from_KT


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

        if chromosome_itr.name not in history_dict:
            path_list.append(Path(Arm(segment_list, "KT_path"), chromosome_itr.name, chromosome_itr.name[:-1]))
            continue

        # label SV segments from history
        history_list = history_dict[chromosome_itr.name]
        for history_entry in history_list:
            event_type = history_entry[0]
            if event_type == "deletion":
                # TODO: add ghost node here
                continue  # cannot label deleted segment

            history_segments = history_entry[1].segments
            for event_segment_itr in history_segments:
                matching_segment_indices = []
                for segment_index in range(len(segment_list)):
                    if (segment_list[segment_index]).same_segment_ignore_dir(event_segment_itr):
                        matching_segment_indices.append(segment_index)

                if len(matching_segment_indices) == 0:
                    raise RuntimeError('history mentioned segment, but segment is not found')
                elif len(matching_segment_indices) == 1:
                    segment_list[matching_segment_indices[0]].segment_type = \
                        "SV" + str(history_counter) + "-" + event_type
                else:
                    segment_noted = False
                    # if more than one segment matched, label the first occurrence that has breakpoint
                    for segment_index in matching_segment_indices:
                        if segment_index == 0:
                            left_cont = True
                        else:
                            left_cont = False
                        if segment_index == len(segment_list) - 1:
                            right_cont = True
                        else:
                            right_cont = False

                        if not left_cont:
                            left_cont = \
                                segments_are_continuous(segment_list[segment_index - 1], segment_list[segment_index])
                        if not right_cont:
                            right_cont = \
                                segments_are_continuous(segment_list[segment_index], segment_list[segment_index + 1])

                        if (not left_cont) or (not right_cont):
                            segment_list[segment_index].segment_type = \
                                "SV" + str(history_counter) + "-" + event_type
                            segment_noted = True
                            break
                    if not segment_noted:
                        raise RuntimeError('SV failed to be labeled')

                history_counter += 1
        path_list.append(Path(Arm(segment_list, "KT_path"), chromosome_itr.name, chromosome_itr.name[:-1]))

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


def test():
    return_list = read_KT_to_path("/media/zhaoyang-new/workspace/KarSim/0908_genomes/KT/23Xe10_r1.kt.txt",
                                  "../Metadata/merged_masking_unique.bed")
    for path in return_list:
        print(path.concise_str())


if __name__ == "__main__":
    test()
