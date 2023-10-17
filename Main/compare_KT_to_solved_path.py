from Structures import *
from read_solved_path import *
from read_KT_to_path import *


def bin_chromosomes(solved_path_file, kt_file, masking_file):
    kt_path_list = read_KT_to_path(kt_file, masking_file)
    solved_path_list = read_solved_path(solved_path_file)
    solved_path_list = rotate_and_bin_path(solved_path_list, masking_file)

    kt_path_bin = {f"Chr{i}": [] for i in range(1, 23)}
    kt_path_bin["ChrX"] = []
    kt_path_bin["ChrY"] = []

    solved_path_bin = {f"Chr{i}": [] for i in range(1, 23)}
    solved_path_bin["ChrX"] = []
    solved_path_bin["ChrY"] = []

    for path in kt_path_list:
        kt_path_bin[path.path_chr].append(path)

    problematic_solve_paths = []
    for path in solved_path_list:
        if path.path_chr.startswith("no centromere") or path.path_chr.startswith("multiple centromere"):
            problematic_solve_paths.append(path)
        else:
            solved_path_bin[path.path_chr].append(path)

    print("\npath with problematic centromere are found (empty if not found): ")
    for path in problematic_solve_paths:
        print(path)

    # check bin sizes
    bin_with_matching_size = []
    bin_with_nonmatching_size = []
    for key in kt_path_bin:
        if len(kt_path_bin[key]) == len(solved_path_bin[key]):
            bin_with_matching_size.append(key)
        else:
            bin_with_nonmatching_size.append(key)

    print("\nbin with non-matching size (empty if not found): ")
    print(bin_with_nonmatching_size)
    print()

    for bin_key in bin_with_matching_size:
        # iterate through all combinations of kt_path vs. solved_path
        # find best pair, pop, and continue finding the next pair
        while len(kt_path_bin[bin_key]) > 0:
            best_score = float('-inf')
            best_kt_alignment = None
            best_solved_path_alignment = None
            best_kt_path_index = None
            best_solved_path_index = None
            for kt_path_index in range(len(kt_path_bin[bin_key])):
                for solved_path_index in range(len(solved_path_bin[bin_key])):
                    current_kt_path = kt_path_bin[bin_key][kt_path_index].duplicate()
                    current_solved_path = solved_path_bin[bin_key][solved_path_index].duplicate()
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
            print("alignment between {}, {}: {}".format(kt_path_bin[bin_key][best_kt_path_index].path_name,
                  solved_path_bin[bin_key][best_solved_path_index].path_name,
                  best_score))
            print(best_kt_alignment)
            print(best_solved_path_alignment)
            print()

            kt_path_bin[bin_key].pop(best_kt_path_index)
            solved_path_bin[bin_key].pop(best_solved_path_index)


def align_paths(segment_list1, segment_list2):
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
            down_value = scoring_matrix[row_index - 1][
                             col_index] + len(segment_list1[row_index - 1]) * indel_penalty_per_nt
            right_value = scoring_matrix[row_index][
                              col_index - 1] + len(segment_list2[col_index - 1]) * indel_penalty_per_nt
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


def test():
    bin_chromosomes(
        "/media/zhaoyang-new/workspace/KarSim/OMKar_outputs/simulation_final/1q21-1_recurrent_microdeletion_r1.1/1q21-1_recurrent_microdeletion_r1.1.txt",
        "/media/zhaoyang-new/workspace/KarSim/1011_genomes/KT/1q21-1_recurrent_microdeletion_r1.kt.txt",
        "Metadata/merged_masking_unique.bed")


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
    test()
