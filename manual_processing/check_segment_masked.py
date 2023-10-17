from Main.Structures import *
from Main.read_masking_regions import read_masking_regions


def check_regions_masked(input_arm: Arm, masking_file):
    masking_arm = read_masking_regions(masking_file)
    if input_arm.arm_intersection(masking_arm):
        return input_arm.report_arm_intersection(masking_arm)
    else:
        return False


def test():
    return_value = check_regions_masked(Arm([Segment('Chr2', 147400000, 148000000)], 'test_arm'),
                                        '../Metadata/merged_masking_unique.bed')
    print(return_value)


def one_time_usage():
    return_value = check_regions_masked(Arm([Segment('Chr15', 74120302, 75680570)], 'test_arm'),
                                        '../Metadata/merged_masking_unique.bed')
    print(return_value)

    return_value = check_regions_masked(Arm([Segment('Chr12', 64678139, 68251745)], 'test_arm'),
                                        '../Metadata/merged_masking_unique.bed')
    print(return_value)

    return_value = check_regions_masked(Arm([Segment('Chr17', 30780079, 31936302)], 'test_arm'),
                                        '../Metadata/merged_masking_unique.bed')
    print(return_value)

    return_value = check_regions_masked(Arm([Segment('Ch17', 36459259, 37856298)], 'test_arm'),
                                        '../Metadata/merged_masking_unique.bed')
    print(return_value)

    return_value = check_regions_masked(Arm([Segment('Chr17', 14194598, 15567589)], 'test_arm'),
                                        '../Metadata/merged_masking_unique.bed')
    print(return_value)

    return_value = check_regions_masked(Arm([Segment('Chr1', 147061832, 148411223)], 'test_arm'),
                                        '../Metadata/merged_masking_unique.bed')
    print(return_value)


if __name__ == "__main__":
    one_time_usage()
