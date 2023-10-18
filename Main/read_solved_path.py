from Structures import *
from read_masking_regions import read_masking_regions


def read_solved_path(file):
    segment_dict = {}
    path_list = []
    with open(file) as fp_read:
        fp_read.readline()

        for line in fp_read:
            line = line.replace("\n", "").split('\t')

            # documenting segments
            if line[0] == "Segment":
                chr_name = str(line[2])
                if chr_name == "23":
                    chr_name = "X"
                elif chr_name == "24":
                    chr_name = "Y"
                chr_name = "Chr" + chr_name
                start = int(line[3].split(".")[0])
                end = int(line[4].split(".")[0])
                segment_dict[int(line[1])] = Segment(chr_name, start, end, "solved_path")
            elif line[0].startswith("Path"):
                line = line[0].split(" = ")
                path_name = line[0]
                line = line[1]
                line = line[:-1].split(" ")
                path_segments = []
                for segment_index_itr in line:
                    direction = segment_index_itr[-1]
                    segment_index_itr = int(segment_index_itr[:-1])
                    new_segment = segment_dict[segment_index_itr].duplicate()
                    if direction == "+":
                        path_segments.append(new_segment)
                    elif direction == "-":
                        new_segment.invert()
                        path_segments.append(new_segment)
                    else:
                        raise ValueError("direction must be + or -")
                path_list.append(Path(Arm(path_segments, "solved_path"), path_name))

    return path_list


def get_centromere_segments(masking_file):
    masking_arm = read_masking_regions(masking_file)
    centromere_segments = []
    for segment in masking_arm.segments:
        if segment.segment_type == "centromere":
            centromere_segments.append(segment)
    return Arm(centromere_segments, "centromeres")


def rotate_and_bin_path(path_list, masking_file):
    """
    only works if each path contains exactly one centromere
    :param masking_file:
    :param path_list:
    :return: path_list
    """
    centromere_arm = get_centromere_segments(masking_file)
    centromere_path = Path(centromere_arm)

    # isolate centromere
    for path in path_list:
        path.generate_mutual_breakpoints(other_path=centromere_path, mutual=False)

    # get centromere, rotate if backward, and bin path
    for path in path_list:
        path_centromere = []
        for centromere_segment in centromere_path.linear_path.segments:
            for segment_itr in path.linear_path.segments:
                if segment_itr.same_segment_ignore_dir(centromere_segment):
                    # if path_centromere is not None:
                    #     raise ValueError("di-centromeric path found")
                    path_centromere.append(segment_itr)
        # if path_centromere is None:
        #     raise ValueError("a-centromeric path found")

        if len(path_centromere) == 1:
            # reverse whole path if centromere found backward
            if not path_centromere[0].direction():
                reversed_segment_list = []
                for segment_itr in reversed(path.linear_path.segments):
                    new_segment = segment_itr.duplicate()
                    new_segment.invert()
                    reversed_segment_list.append(new_segment)
                path.linear_path.segments = reversed_segment_list
            # assign bin
            path.path_chr = path_centromere[0].chr_name
        elif len(path_centromere) == 0:
            path.path_chr = "no centromere"
        else:
            path.path_chr = "multiple centromere: "
            for centromere_itr in path_centromere:
                path.path_chr += centromere_itr.chr_name + " "

    return path_list


def report_centromere_anomaly(path_list):
    print(len(path_list))
    for path in path_list:
        if path.path_chr.startswith("no centromere") or path.path_chr.startswith("multiple centromere"):
            print(path.get_path_notes())


def test():
    path_list = read_solved_path("/media/zhaoyang-new/workspace/KarSim/OMKar_outputs/simulation_final/1q21-1_recurrent_microdeletion_r1.1/1q21-1_recurrent_microdeletion_r1.1.txt")
    path_list = rotate_and_bin_path(path_list, "../Metadata/merged_masking_unique.bed")
    # report_centromere_anomaly(path_list)
    for path in path_list:
        print(path)


def cmd():
    import argparse
    parser = argparse.ArgumentParser(description="dicentromeric and acentromeric checker")
    parser.add_argument("--file", type=str, dest='omkar_file', help="file path to OMKar's solved path")
    args = parser.parse_args()

    path_list = read_solved_path(args.omkar_file)
    for path in path_list:
        path.linear_path.merge_breakpoints()
    path_list = rotate_and_bin_path(path_list, "Metadata/merged_masking_unique.bed")
    report_centromere_anomaly(path_list)


if __name__ == "__main__":
    cmd()
