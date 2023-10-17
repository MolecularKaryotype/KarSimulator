from Structures import *
from read_masking_regions import read_masking_regions
from Start_Genome import generate_genome_from_KT


def read_KT_to_path(KT_file, masking_file):
    genome = generate_genome_from_KT(KT_file)
    t2_segments = get_t2_segments(masking_file)
    path_list = []

    for chromosome_itr in genome:
        t1_segment = Segment(chromosome_itr.name[:-1], 0, chromosome_itr.t1_len - 1, "telomere1")
        segment_list = [t1_segment]

        for segment_itr in chromosome_itr.p_arm.segments:
            segment_itr.segment_type = "p_arm"
            segment_list.append(segment_itr)

        # TODO: look into why centromere segment is not an instance of Segment: isinstance() = False
        cen_segment = chromosome_itr.centromere.segments[0].duplicate()
        cen_segment.segment_type = "centromere"
        segment_list.append(cen_segment)

        for segment_itr in chromosome_itr.q_arm.segments:
            segment_itr.segment_type = "q_arm"
            segment_list.append(segment_itr)

        t2_segment = None
        for find_t2_segment_itr in t2_segments.segments:
            if find_t2_segment_itr.chr_name == chromosome_itr.name[:-1]:
                t2_segment = find_t2_segment_itr
                break
        segment_list.append(t2_segment)

        path_list.append(Path(Arm(segment_list, "KT_path"), chromosome_itr.name, chromosome_itr.name[:-1]))

    return path_list


def get_t2_segments(masking_file):
    masking_arm = read_masking_regions(masking_file)
    t2_segments = []
    for segment in masking_arm.segments:
        if segment.segment_type == "telomere2":
            t2_segments.append(segment)
    return Arm(t2_segments, "t2_segments")


def test():
    return_list = read_KT_to_path("/media/zhaoyang-new/workspace/KarSim/1011_genomes/KT/1q21-1_recurrent_microdeletion_r1.kt.txt",
                                  "../Metadata/merged_masking_unique.bed")
    for path in return_list:
        print(path)


if __name__ == "__main__":
    test()
