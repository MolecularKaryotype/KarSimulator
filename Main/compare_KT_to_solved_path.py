from Structures import Segment, Arm, Path
from read_solved_path import *
from read_KT_to_path import *
from read_masking_regions import read_masking_regions
from Start_Genome import generate_genome_from_KT


def compare_paths(solved_path_file, kt_file, masking_file):
    kt_path_list = read_KT_to_path(kt_file, masking_file)
    genome = generate_genome_from_KT(kt_file, ordinal_info_included=True)
    solved_path_list = read_solved_path(solved_path_file)

    label_acrocentric_regions('../Metadata/comparison_forbidden_regions.bed', solved_path_list)
    label_acrocentric_regions('../Metadata/comparison_forbidden_regions.bed', kt_path_list)

    # create bins of dependent chrs
    bin_list = bin_dependent_chr(kt_path_list)

    # debugging a specific alignment
    # current_kt_path = kt_path_list[16].duplicate()
    # current_solved_path = solved_path_list[12].duplicate()
    # current_kt_path.generate_mutual_breakpoints(other_path=current_solved_path, mutual=True)
    # best_score, best_kt_alignment, best_solved_path_alignment, sv_total, sv_captured, jaccard = \
    #     align_paths(current_kt_path.linear_path.segments, current_solved_path.linear_path.segments)
    # print("best score: " + str(best_score))
    # print(best_kt_alignment)
    # print(best_solved_path_alignment)

    chr_pairing = {}
    while len(kt_path_list) > len(chr_pairing) and len(solved_path_list) > len(chr_pairing.values()):
        best_score = float('-inf')
        best_kt_index = None
        best_path_index = None
        best_path_is_reversed = False
        for kt_path_index in range(len(kt_path_list)):
            if kt_path_index in chr_pairing:
                continue
            for solved_path_index in range(len(solved_path_list)):
                if solved_path_index in chr_pairing.values():
                    continue
                current_kt_path = kt_path_list[kt_path_index].duplicate()
                current_solved_path = solved_path_list[solved_path_index].duplicate()
                current_kt_path.generate_mutual_breakpoints(other_path=current_solved_path, mutual=True)

                current_score, current_kt_alignment, current_solved_path_alignment, sv_total, sv_captured, jaccard = \
                    align_paths(current_kt_path.linear_path.segments,
                                current_solved_path.linear_path.segments)
                if current_score > best_score:
                    best_score = current_score
                    best_kt_index = kt_path_index
                    best_path_index = solved_path_index
                    best_path_is_reversed = False

                # check if reversed gives higher score
                current_kt_path = kt_path_list[kt_path_index].duplicate()
                current_solved_path = solved_path_list[solved_path_index].duplicate()
                current_solved_path.reverse()
                current_kt_path.generate_mutual_breakpoints(other_path=current_solved_path, mutual=True)

                current_score, current_kt_alignment, current_solved_path_alignment, sv_total, sv_captured, jaccard = \
                    align_paths(current_kt_path.linear_path.segments,
                                current_solved_path.linear_path.segments)
                if current_score > best_score:
                    best_score = current_score
                    best_kt_index = kt_path_index
                    best_path_index = solved_path_index
                    best_path_is_reversed = True

        # record best pairing
        if best_path_is_reversed:
            solved_path_list[best_path_index].reverse()
        chr_pairing[best_kt_index] = best_path_index

    # summary statistics
    print("unaligned standard: (empty if none)")
    if len(kt_path_list) > len(chr_pairing):
        for path_index in range(len(kt_path_list)):
            if path_index not in chr_pairing:
                print(kt_path_list[path_index])
    print("unaligned reconstruction: (empty if none)")
    if len(solved_path_list) > len(chr_pairing.values()):
        for path_index in range(len(solved_path_list)):
            if path_index not in chr_pairing.values():
                print(solved_path_list[path_index])

    genome_output_str = ""
    bin_jaccard_strs = []
    genome_sv_captured = {}
    genome_sv_total = {}
    genome_indel = 0
    for bin_index in range(len(bin_list)):
        bin_output_str = ""
        bin_sv_captured = {}
        bin_sv_total = {}
        bin_indel = 0

        chr_list_str = ",".join(bin_list[bin_index][0])
        genome_output_str += "\ndependent component " + str(bin_index) + ": " + chr_list_str + "\n"
        for kt_path_index_itr in bin_list[bin_index][1]:
            # record output for current pair
            if kt_path_index_itr in chr_pairing:
                current_solved_path_index = chr_pairing[kt_path_index_itr]
            else:
                bin_output_str += "kt unaligned: " + str(kt_path_list[kt_path_index_itr]) + "\n"
                continue

            # re-perform alignment to get statistics
            current_kt_path = kt_path_list[kt_path_index_itr].duplicate()
            current_solved_path = solved_path_list[current_solved_path_index].duplicate()
            current_kt_path.generate_mutual_breakpoints(other_path=current_solved_path, mutual=True)
            best_score, best_kt_alignment, best_solved_path_alignment, sv_total, sv_captured, jaccard = \
                align_paths(current_kt_path.linear_path.segments,
                            current_solved_path.linear_path.segments)

            # append alignment for output
            bin_output_str += "alignment between {}, {}: {}\n".format(kt_path_list[kt_path_index_itr].path_name,
                                                                      solved_path_list[current_solved_path_index].path_name,
                                                                      best_score)
            bin_output_str += "Alignment's Jaccard score: {}\n".format(jaccard)
            bin_output_str += best_kt_alignment + "\n"
            bin_output_str += best_solved_path_alignment + "\n"

            # append stats to dependent component (bin) level
            bin_indel += best_score
            for sv_name in sv_total:
                if sv_name not in bin_sv_total:
                    bin_sv_total[sv_name] = 0
                    bin_sv_captured[sv_name] = 0
                bin_sv_total[sv_name] += sv_total[sv_name]
                bin_sv_captured[sv_name] += sv_captured[sv_name]
                # DEBUG
                # print("sv_captured[sv_name]: " + str(sv_captured[sv_name]))
                # print("bin_sv_captured[sv_name]: " + str(bin_sv_captured[sv_name]))
                # print(bin_sv_total[sv_name])

        bin_sv_str = ""
        bin_jaccard_num = 0
        bin_jaccard_denom = 0
        for sv_name in bin_sv_total:
            sv_capture_rate = bin_sv_captured[sv_name] / bin_sv_total[sv_name]
            sv_capture_rate = "{:.4f}".format(sv_capture_rate)
            sv_history_info = genome.history[int(sv_name.replace('SV', ''))]
            sv_history_str = '{' + sv_history_info[0] + ' from ' + sv_history_info[2].name + " to " + sv_history_info[3].name + '}'
            bin_sv_str += "SV_name: {}-{}\tSV_total: {}\tSV_captured: {}\tSV_capture_rate: {}\n".format(sv_name,
                                                                                                            sv_history_str,
                                                                                                            str(bin_sv_total[sv_name]),
                                                                                                            str(bin_sv_captured[sv_name]),
                                                                                                            sv_capture_rate)
            bin_jaccard_num += bin_sv_captured[sv_name]
            bin_jaccard_denom += bin_sv_total[sv_name]
        bin_jaccard_denom += abs(bin_indel)
        bin_jaccard = bin_jaccard_num / bin_jaccard_denom
        bin_jaccard_str = "component {} ({}), #SV: {}, Jaccard Score: ".format(str(bin_index), chr_list_str, len(bin_sv_total))
        if len(bin_sv_total) == 0:
            bin_jaccard_str += 'No Event\n'
        else:
            bin_jaccard_str += str(bin_jaccard) + '\n'
        genome_output_str += bin_jaccard_str + bin_sv_str + "\n" + bin_output_str + "\n"
        if len(bin_sv_total) != 0:
            bin_jaccard_strs.append(bin_jaccard_str)

        # append stats to genome level
        genome_indel += bin_indel
        for sv_name in bin_sv_total:
            if sv_name not in genome_sv_total:
                genome_sv_total[sv_name] = 0
                genome_sv_captured[sv_name] = 0
            genome_sv_total[sv_name] += bin_sv_total[sv_name]
            genome_sv_captured[sv_name] += bin_sv_captured[sv_name]

    def custom_sort_sv_name(key):
        return int(key[2:])
    sorted_keys = sorted(genome_sv_total.keys(), key=custom_sort_sv_name)
    genome_sv_str = ""
    sv_captured_number = 0
    genome_jaccard_num = 0
    genome_jaccard_denom = 0
    for sv_name in sorted_keys:
        sv_capture_rate = genome_sv_captured[sv_name] / genome_sv_total[sv_name]
        if sv_capture_rate > 0.5:
            sv_captured_number += 1
        sv_capture_rate = "{:.4f}".format(sv_capture_rate)
        sv_history_info = genome.history[int(sv_name.replace('SV', ''))]
        sv_history_str = '{' + sv_history_info[0] + ' from ' + sv_history_info[2].name + " to " + sv_history_info[3].name + '}'
        genome_sv_str += "SV_name: {}-{}\tSV_total: {}\tSV_captured: {}\tSV_capture_rate: {}\n".format(sv_name,
                                                                                                       sv_history_str,
                                                                                                       str(genome_sv_total[sv_name]),
                                                                                                       str(genome_sv_captured[sv_name]),
                                                                                                       sv_capture_rate)
        genome_jaccard_num += genome_sv_captured[sv_name]
        genome_jaccard_denom += genome_sv_total[sv_name]

    genome_jaccard_denom += abs(genome_indel)
    genome_jaccard = genome_jaccard_num / genome_jaccard_denom
    print("genome indel: {}".format(abs(genome_indel)))
    print("genome Jaccard Score: {}".format(str(genome_jaccard)))
    print("connected components' Jaccard Score: ")
    print(''.join(bin_jaccard_strs))
    print("SV: {} present, {} captured".format(str(len(sorted_keys)), str(sv_captured_number)))
    print(genome_sv_str)
    print(genome_output_str)


def bin_dependent_chr(kt_path_list):
    class Bin:
        def __init__(self):
            self.chr_components = set()
            self.kt_indices = []

        def chr_in(self, input_chr):
            return input_chr in self.chr_components

        def chr_add(self, input_chr):
            self.chr_components.add(input_chr)

        def index_append(self, input_index):
            self.kt_indices.append(input_index)

        def output_indices(self):
            return self.kt_indices

    bin_list = []
    for kt_path_index in range(len(kt_path_list)):
        current_chr = set()
        for segment in kt_path_list[kt_path_index].linear_path.segments:
            current_chr.add(segment.chr_name)

        bin_found = None
        for bin_itr in bin_list:
            for chr_itr in current_chr:
                if bin_itr.chr_in(chr_itr):
                    bin_found = bin_itr
                    break
        if bin_found is None:
            new_bin = Bin()
            new_bin.chr_components = current_chr
            new_bin.index_append(kt_path_index)
            bin_list.append(new_bin)
        else:
            for chr_itr in current_chr:
                bin_found.chr_add(chr_itr)
            bin_found.index_append(kt_path_index)

    return_list = []
    for bin_itr in bin_list:
        return_list.append(tuple([bin_itr.chr_components, bin_itr.output_indices()]))
    return return_list


def label_acrocentric_regions(acrocentric_forbidden_file, input_paths):
    acrocentric_region = read_masking_regions(acrocentric_forbidden_file)
    for path in input_paths:
        tmp_acrocentric_path = Path(acrocentric_region.duplicate(), "forbidden_regions")
        path.generate_mutual_breakpoints(other_path=tmp_acrocentric_path, mutual=True)
        for segment in path.linear_path.segments:
            for acrocentric_segment_itr in tmp_acrocentric_path.linear_path.segments:
                if segment.same_segment_ignore_dir(acrocentric_segment_itr):
                    segment.segment_type = "acrocentric"
                    break


def label_telomere_regions(all_forbidden_region_file, input_paths):
    all_forbidden_region = read_masking_regions(all_forbidden_region_file)
    telomere_segments = []
    for segment in all_forbidden_region.segments:
        if segment.segment_type in ['telomere1', 'telomere2']:
            telomere_segments.append(segment)
    telomere_regions = Arm(telomere_segments, "telomeres")

    for path in input_paths:
        tmp_telomere_path = Path(telomere_regions.duplicate(), "forbidden_regions")
        path.generate_mutual_breakpoints(other_path=tmp_telomere_path, mutual=True)
        for segment in path.linear_path.segments:
            for telomere_segment_itr in tmp_telomere_path.linear_path.segments:
                if segment.same_segment_ignore_dir(telomere_segment_itr):
                    segment.segment_type = telomere_segment_itr.segment_type
                    break


def align_paths(segment_list1, segment_list2):
    forbidden_comparison_region_types = ['acrocentric', 'telomere1', 'telomere2']
    indel_penalty_per_nt = -1
    alignment_1 = []
    alignment_2 = []
    scoring_matrix = [[0 for i in range(0, len(segment_list2) + 1)] for j in
                      range(0, len(segment_list1) + 1)]
    backtrack_matrix = [["" for i in range(0, len(segment_list2) + 1)] for j in
                        range(0, len(segment_list1) + 1)]

    # initialize starting grid
    scoring_matrix[0][0] = 0
    for row_index in range(1, len(segment_list1) + 1):
        current_segment = segment_list1[row_index - 1]
        if current_segment.segment_type is None or current_segment.segment_type not in ['telomere1', 'telomere2', 'centromere', 'acrocentric']:
            scoring_matrix[row_index][0] = scoring_matrix[row_index - 1][0] + len(current_segment) * indel_penalty_per_nt
        else:
            scoring_matrix[row_index][0] = scoring_matrix[row_index - 1][0]
        backtrack_matrix[row_index][0] = "down"
    for col_index in range(1, len(segment_list2) + 1):
        current_segment = segment_list2[col_index - 1]
        if current_segment.segment_type is None or current_segment.segment_type not in ['telomere1', 'telomere2', 'centromere', 'acrocentric']:
            scoring_matrix[0][col_index] = scoring_matrix[0][col_index - 1] + len(current_segment) * indel_penalty_per_nt
        else:
            scoring_matrix[0][col_index] = scoring_matrix[0][col_index - 1]
        backtrack_matrix[0][col_index] = "rigt"

    for row_index in range(1, len(segment_list1) + 1):
        for col_index in range(1, len(segment_list2) + 1):
            if segment_list1[row_index - 1].segment_type is not None and segment_list1[row_index - 1].segment_type in forbidden_comparison_region_types:
                # forbidden region indel
                down_value = scoring_matrix[row_index - 1][col_index]
            elif segment_list1[row_index - 1].segment_type is not None and segment_list1[row_index - 1].segment_type.startswith('del'):
                # ghost indel
                down_value = scoring_matrix[row_index - 1][col_index]
            else:
                # standard indel
                down_value = scoring_matrix[row_index - 1][col_index] \
                             + len(segment_list1[row_index - 1]) * indel_penalty_per_nt

            if segment_list2[col_index - 1].segment_type is not None and segment_list2[col_index - 1].segment_type in forbidden_comparison_region_types:
                # forbidden region indel
                right_value = scoring_matrix[row_index][col_index - 1]
            else:
                # standard indel
                right_value = scoring_matrix[row_index][col_index - 1] \
                              + len(segment_list2[col_index - 1]) * indel_penalty_per_nt

            if segment_list1[row_index - 1].segment_type is not None and segment_list1[row_index - 1].segment_type.startswith('del'):
                # matching a ghost segment
                diagonal_value = scoring_matrix[row_index - 1][col_index - 1] + len(segment_list1[row_index - 1]) * indel_penalty_per_nt
            elif segment_list1[row_index - 1] == segment_list2[col_index - 1]:
                # match
                diagonal_value = scoring_matrix[row_index - 1][col_index - 1]
            else:
                # mismatch: not allowed
                diagonal_value = float('-inf')

            if diagonal_value >= down_value and diagonal_value >= right_value:
                scoring_matrix[row_index][col_index] = diagonal_value
                backtrack_matrix[row_index][col_index] = "diag"
            elif down_value >= right_value:
                scoring_matrix[row_index][col_index] = down_value
                backtrack_matrix[row_index][col_index] = "down"
            else:
                scoring_matrix[row_index][col_index] = right_value
                backtrack_matrix[row_index][col_index] = "rigt"

    # backtracking
    final_score = scoring_matrix[len(segment_list1)][len(segment_list2)]
    current_row = len(segment_list1)
    current_col = len(segment_list2)
    sv_total = {}
    sv_captured = {}
    for segment_itr in segment_list1:
        if segment_itr.segment_type is not None and segment_itr.segment_type not in ['telomere1', 'telomere2', 'centromere', 'acrocentric']:
            event_name = ': '.join(segment_itr.segment_type.split(': ')[1:])
            if event_name not in sv_total:
                sv_total[event_name] = 0
            if len(event_name.split(': ')) == 2 and event_name.split(': ')[0] not in sv_total:
                sv_total[event_name.split(': ')[0]] = 0
            sv_total[event_name] += len(segment_itr)
    for event_name in sv_total:
        sv_captured[event_name] = 0

    while True:
        if current_row == 0 and current_col == 0:
            break
        if backtrack_matrix[current_row][current_col] == "diag":
            alignment_1.insert(0, segment_list1[current_row - 1])
            alignment_2.insert(0, segment_list2[current_col - 1])
            if segment_list1[current_row - 1].segment_type is not None and segment_list1[current_row - 1].segment_type not in ['telomere1', 'telomere2',
                                                                                                                               'centromere', 'acrocentric']:
                if not segment_list1[current_row - 1].segment_type.startswith('del'):
                    event_name = ': '.join(segment_list1[current_row - 1].segment_type.split(': ')[1:])
                    sv_captured[event_name] += len(segment_list1[current_row - 1])
            current_col -= 1
            current_row -= 1
        elif backtrack_matrix[current_row][current_col] == "down":
            alignment_1.insert(0, segment_list1[current_row - 1])
            alignment_2.insert(0, "-")
            if segment_list1[current_row - 1].segment_type is not None and segment_list1[current_row - 1].segment_type not in ['telomere1', 'telomere2',
                                                                                                                               'centromere', 'acrocentric']:
                if segment_list1[current_row - 1].segment_type.startswith('del'):
                    event_name = ': '.join(segment_list1[current_row - 1].segment_type.split(': ')[1:])
                    sv_captured[event_name] += len(segment_list1[current_row - 1])
            current_row -= 1
        elif backtrack_matrix[current_row][current_col] == "rigt":
            alignment_1.insert(0, "-")
            alignment_2.insert(0, segment_list2[current_col - 1])
            current_col -= 1
        else:
            print("error in backtrack matrix")

    alignment_1_string = ""
    alignment_2_string = ""
    for index in alignment_1:
        if index == "-":
            alignment_1_string += index
        else:
            alignment_1_string += index.alignment_output()
        alignment_1_string += "\t"
    for index in alignment_2:
        if index == "-":
            alignment_2_string += index
        else:
            alignment_2_string += index.alignment_output()
        alignment_2_string += "\t"

    # remove duplicate ghost segment's contribution in SV total and SV captured
    for sv_name in sv_total:
        if len(sv_name.split(': ')) == 2:
            # it is a deletion with ghost multiplicity
            ghost_multiplicity = int(sv_name.split(': ')[1])
            individual_count = sv_total[sv_name] / ghost_multiplicity
            sv_total[sv_name] = sv_total[sv_name] - (ghost_multiplicity - 1) * individual_count
            sv_captured[sv_name] = sv_captured[sv_name] - (ghost_multiplicity - 1) * individual_count
    # merge ins and del of the same SV
    for sv_name in sv_total:
        if len(sv_name.split(': ')) == 2:
            sv_total[sv_name.split(': ')[0]] += sv_total[sv_name]
            sv_captured[sv_name.split(': ')[0]] += sv_captured[sv_name]
    new_sv_total = {}
    new_sv_captured = {}
    for sv_name in sv_total:
        if len(sv_name.split(': ')) == 1:
            new_sv_total[sv_name] = sv_total[sv_name]
            new_sv_captured[sv_name] = sv_captured[sv_name]
    sv_total = new_sv_total
    sv_captured = new_sv_captured

    # calculate jaccard
    numerator = 0
    denominator = abs(final_score)
    for sv_name in sv_total:
        numerator += sv_captured[sv_name]
        denominator += sv_total[sv_name]
    jaccard = numerator / denominator
    if len(sv_total) == 0:
        jaccard = "No Event"
    else:
        jaccard = "{:.4f}".format(jaccard)
    return final_score, alignment_1_string, alignment_2_string, sv_total, sv_captured, jaccard


def test_compare_paths():
    omkar_file = "../scoring_files/modified_OMKar/23Y_WAGR_11p13_deletion_r2.1.txt"
    kt_file = "../scoring_files/modified_KT/23Y_WAGR_11p13_deletion_r2.kt.txt"
    masking_file = "../Metadata/merged_masking_unique.bed"
    compare_paths(omkar_file, kt_file, masking_file)
    # compare_paths(
    #     "/media/zhaoyang-new/workspace/KarSim/KarSimulator/test_folder/23Xe10_r1.1.txt",
    #     "/media/zhaoyang-new/workspace/KarSim/KarSimulator/test_folder/23Xe10_r1.kt.txt",
    #     "../Metadata/merged_masking_unique.bed")
    # compare_paths(
    #     "/media/zhaoyang-new/workspace/KarSim/OMKar_outputs/simulation_final_v3/23X_22q11-2_distal_deletion_r1.1/23X_22q11-2_distal_deletion_r1.1.txt",
    #     "/media/zhaoyang-new/workspace/KarSim/KarSimulator/scoring_files/23X_22q11-2_distal_deletion_r1.kt.txt",
    #     "../Metadata/merged_masking_unique.bed")
    # compare_paths(
    #     "/media/zhaoyang-new/workspace/KarSim/OMKar_outputs/simulation_final_v3/STS_r1.1/STS_r1.1.txt",
    #     "/media/zhaoyang-new/workspace/KarSim/KarSimulator/scoring_files/STS_r1.kt.txt",
    #     "../Metadata/merged_masking_unique.bed")


def test_align_paths():
    segment_list1 = [Segment("Chr1", 0, 9999, "telomere1"),
                     Segment("Chr1", 10000, 19999, "p_arm"),
                     Segment("Chr1", 20000, 24999, "p_arm"),
                     Segment("Chr1", 25000, 29999, "p_arm"),
                     Segment("Chr1", 30000, 40000, "p_arm")]
    segment_list2 = [Segment("Chr1", 0, 9999, "telomere1"),
                     Segment("Chr1", 10000, 19999, "p_arm"),
                     Segment("Chr1", 20000, 24999, "p_arm"),
                     Segment("Chr1", 30000, 40000, "p_arm")]

    a, b, c = align_paths(segment_list1, segment_list2)
    print(a)
    print(b)
    print(c)


if __name__ == "__main__":
    test_compare_paths()
