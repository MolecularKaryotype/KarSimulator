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

    while len(kt_path_list) > 0 and len(solved_path_list) > 0:
        best_score = float('-inf')
        best_kt_alignment = None
        best_solved_path_alignment = None
        best_kt_path_index = None
        best_solved_path_index = None
        best_sv_total = None
        best_sv_captured = None
        best_jaccard = None
        for kt_path_index in range(len(kt_path_list)):
            for solved_path_index in range(len(solved_path_list)):
                current_kt_path = kt_path_list[kt_path_index].duplicate()
                current_solved_path = solved_path_list[solved_path_index].duplicate()
                current_kt_path.generate_mutual_breakpoints(other_path=current_solved_path, mutual=True)

                current_score, current_kt_alignment, current_solved_path_alignment, sv_total, sv_captured, jaccard = \
                    align_paths(current_kt_path.linear_path.segments,
                                current_solved_path.linear_path.segments)
                if current_score > best_score:
                    best_score = current_score
                    best_kt_alignment = current_kt_alignment
                    best_solved_path_alignment = current_solved_path_alignment
                    best_kt_path_index = kt_path_index
                    best_solved_path_index = solved_path_index
                    best_sv_total = sv_total
                    best_sv_captured = sv_captured
                    best_jaccard = jaccard

                # check if reversed gives higher score
                current_kt_path = kt_path_list[kt_path_index].duplicate()
                current_solved_path = solved_path_list[solved_path_index].duplicate()
                current_solved_path.reverse()
                current_kt_path.generate_mutual_breakpoints(other_path=current_solved_path, mutual=True)

                current_score, current_kt_alignment, current_solved_path_alignment, sv_total, sv_captured, jaccard= \
                    align_paths(current_kt_path.linear_path.segments,
                                current_solved_path.linear_path.segments)
                if current_score > best_score:
                    best_score = current_score
                    best_kt_alignment = current_kt_alignment
                    best_solved_path_alignment = current_solved_path_alignment
                    best_kt_path_index = kt_path_index
                    best_solved_path_index = solved_path_index
                    best_sv_total = sv_total
                    best_sv_captured = sv_captured
                    best_jaccard = jaccard

        # output current pair
        print("alignment between {}, {}: {}".format(kt_path_list[best_kt_path_index].path_name,
              solved_path_list[best_solved_path_index].path_name, best_score))
        for sv_name in best_sv_total:
            sv_capture_rate = best_sv_captured[sv_name] / best_sv_total[sv_name]
            sv_capture_rate = "{:.4f}".format(sv_capture_rate)
            sv_history_info = genome.history[int(sv_name.replace('SV', ''))]
            sv_history_str = '{' + sv_history_info[0] + ' from ' + sv_history_info[2].name + " to " + sv_history_info[3].name + '}'
            print("SV_name: {}-{}\tSV_total: {}\tSV_captured: {}\tSV_capture_rate: {}".format(sv_name, sv_history_str, str(best_sv_total[sv_name]), str(best_sv_captured[sv_name]), sv_capture_rate))
        print("jaccard score: {}".format(best_jaccard))
        print(best_kt_alignment)
        print(best_solved_path_alignment)
        print()

        kt_path_list.pop(best_kt_path_index)
        solved_path_list.pop(best_solved_path_index)

    if len(kt_path_list) > 0:
        for path in kt_path_list:
            print(path)
    if len(solved_path_list) > 0:
        for path in solved_path_list:
            print(path)


def label_acrocentric_regions(acrocentric_forbidden_file, input_paths):
    acrocentric_region = read_masking_regions(acrocentric_forbidden_file)
    for path in input_paths:
        tmp_acrocentric_path = Path(acrocentric_region.duplicate(), "forbidden_regions")
        path.generate_mutual_breakpoints(other_path=tmp_acrocentric_path, mutual=True)
        for segment in path.linear_path.segments:
            if segment in tmp_acrocentric_path.linear_path.segments:
                segment.segment_type = "acrocentric"


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
        scoring_matrix[row_index][0] = 0
        backtrack_matrix[row_index][0] = "down"
    for col_index in range(1, len(segment_list2) + 1):
        scoring_matrix[0][col_index] = 0
        backtrack_matrix[0][col_index] = "rigt"

    for row_index in range(1, len(segment_list1) + 1):
        for col_index in range(1, len(segment_list2) + 1):
            if segment_list1[row_index - 1].segment_type in forbidden_comparison_region_types:
                down_value = scoring_matrix[row_index - 1][col_index]
            elif segment_list1[row_index - 1].segment_type is not None and segment_list1[row_index - 1].segment_type.startswith('del'):
                down_value = scoring_matrix[row_index - 1][col_index]
            else:
                down_value = scoring_matrix[row_index - 1][col_index] \
                             + len(segment_list1[row_index - 1]) * indel_penalty_per_nt

            if segment_list2[col_index - 1].segment_type in forbidden_comparison_region_types:
                right_value = scoring_matrix[row_index][col_index - 1]
            else:
                right_value = scoring_matrix[row_index][col_index - 1] \
                              + len(segment_list2[col_index - 1]) * indel_penalty_per_nt

            if segment_list1[row_index - 1] == segment_list2[col_index - 1]:
                diagonal_value = scoring_matrix[row_index - 1][col_index - 1]
            elif segment_list1[row_index - 1].segment_type is not None and segment_list1[row_index - 1].segment_type.startswith('del'):
                diagonal_value = scoring_matrix[row_index - 1][col_index - 1] + + len(segment_list1[row_index - 1]) * indel_penalty_per_nt
            else:
                diagonal_value = float('-inf')  # mismatch not allowed

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
            if segment_list1[current_row - 1].segment_type is not None and segment_list1[current_row - 1].segment_type not in ['telomere1', 'telomere2', 'centromere', 'acrocentric']:
                if not segment_list1[current_row - 1].segment_type.startswith('del'):
                    event_name = ': '.join(segment_list1[current_row - 1].segment_type.split(': ')[1:])
                    sv_captured[event_name] += len(segment_list1[current_row - 1])
            current_col -= 1
            current_row -= 1
        elif backtrack_matrix[current_row][current_col] == "down":
            alignment_1.insert(0, segment_list1[current_row - 1])
            alignment_2.insert(0, "-")
            if segment_list1[current_row - 1].segment_type is not None and segment_list1[current_row - 1].segment_type not in ['telomere1', 'telomere2', 'centromere', 'acrocentric']:
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
            alignment_1_string += str(index)
        alignment_1_string += "\t"
    for index in alignment_2:
        if index == "-":
            alignment_2_string += index
        else:
            alignment_2_string += str(index)
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
    jaccard = "{:.4f}".format(jaccard)
    return final_score, alignment_1_string, alignment_2_string, sv_total, sv_captured, jaccard


def test_compare_paths():
    compare_paths(
        "/media/zhaoyang-new/workspace/KarSim/OMKar_outputs/simulation_final_v3/23Y_Cri_du_Chat_r1.1/23Y_Cri_du_Chat_r1.1.txt",
        "/media/zhaoyang-new/workspace/KarSim/KarSimulator/modified_KT/23Y_Cri_du_Chat_r1.kt.txt",
        "../Metadata/merged_masking_unique.bed")
    # compare_paths(
    #     "/media/zhaoyang-new/workspace/KarSim/KarSimulator/test_folder/23Xe10_r1.1.txt",
    #     "/media/zhaoyang-new/workspace/KarSim/KarSimulator/test_folder/23Xe10_r1.kt.txt",
    #     "../Metadata/merged_masking_unique.bed")
    # compare_paths(
    #     "/media/zhaoyang-new/workspace/KarSim/OMKar_outputs/simulation_final_v3/23X_22q11-2_distal_deletion_r1.1/23X_22q11-2_distal_deletion_r1.1.txt",
    #     "/media/zhaoyang-new/workspace/KarSim/KarSimulator/modified_KT/23X_22q11-2_distal_deletion_r1.kt.txt",
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
