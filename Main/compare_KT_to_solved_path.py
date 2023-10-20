from Structures import Segment, Arm, Path
from read_solved_path import *
from read_KT_to_path import *
from read_masking_regions import read_masking_regions


def compare_paths(solved_path_file, kt_file, masking_file):
    kt_path_list = read_KT_to_path(kt_file, masking_file)
    solved_path_list = read_solved_path(solved_path_file)

    label_acrocentric_regions('../Metadata/comparison_forbidden_regions.bed', solved_path_list)
    label_acrocentric_regions('../Metadata/comparison_forbidden_regions.bed', kt_path_list)

    while len(kt_path_list) > 0 and len(solved_path_list) > 0:
        best_score = float('-inf')
        best_kt_alignment = None
        best_solved_path_alignment = None
        best_kt_path_index = None
        best_solved_path_index = None
        for kt_path_index in range(len(kt_path_list)):
            for solved_path_index in range(len(solved_path_list)):
                current_kt_path = kt_path_list[kt_path_index].duplicate()
                current_solved_path = solved_path_list[solved_path_index].duplicate()
                current_kt_path.generate_mutual_breakpoints(other_path=current_solved_path, mutual=True)

                current_score, current_kt_alignment, current_solved_path_alignment = \
                    align_paths(current_kt_path.linear_path.segments,
                                current_solved_path.linear_path.segments)
                if current_score > best_score:
                    best_score = current_score
                    best_kt_alignment = current_kt_alignment
                    best_solved_path_alignment = current_solved_path_alignment
                    best_kt_path_index = kt_path_index
                    best_solved_path_index = solved_path_index

                # check if reversed gives higher score
                current_kt_path = kt_path_list[kt_path_index].duplicate()
                current_solved_path = solved_path_list[solved_path_index].duplicate()
                current_solved_path.reverse()
                current_kt_path.generate_mutual_breakpoints(other_path=current_solved_path, mutual=True)

                current_score, current_kt_alignment, current_solved_path_alignment = \
                    align_paths(current_kt_path.linear_path.segments,
                                current_solved_path.linear_path.segments)
                if current_score > best_score:
                    best_score = current_score
                    best_kt_alignment = current_kt_alignment
                    best_solved_path_alignment = current_solved_path_alignment
                    best_kt_path_index = kt_path_index
                    best_solved_path_index = solved_path_index

        # output current pair
        print("alignment between {}, {}: {}".format(kt_path_list[best_kt_path_index].path_name,
              solved_path_list[best_solved_path_index].path_name, best_score))
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

    while True:
        if current_row == 0 and current_col == 0:
            break
        if backtrack_matrix[current_row][current_col] == "diag":
            alignment_1.insert(0, segment_list1[current_row - 1])
            alignment_2.insert(0, segment_list2[current_col - 1])
            current_col -= 1
            current_row -= 1
        elif backtrack_matrix[current_row][current_col] == "down":
            alignment_1.insert(0, segment_list1[current_row - 1])
            alignment_2.insert(0, "-")
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
    return final_score, alignment_1_string, alignment_2_string


def test_compare_paths():
    compare_paths(
        "/Users/zhaoyangjia/Library/CloudStorage/OneDrive-UCSanDiego/bionano/Solved_paths/simulation_final/1q21-1_recurrent_microdeletion_r1.1/1q21-1_recurrent_microdeletion_r1.1.txt",
        "/Users/zhaoyangjia/Library/CloudStorage/OneDrive-UCSanDiego/bionano/Solved_paths/simulation_final/1q21-1_recurrent_microdeletion_r1.1/1q21-1_recurrent_microdeletion_r1.kt.txt",
        "../Metadata/merged_masking_unique.bed")



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
