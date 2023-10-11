from Main.Structures import *
from Main.read_regions_regions import read_masking_regions


def check_regions_masked(input_arm: Arm, masking_file):
    masking_arm = read_masking_regions(masking_file)
    return input_arm.arm_intersection(masking_arm)


def test():
    return_value = check_regions_masked(Arm([Segment('Chr2', 147400000, 148000000)], 'test_arm'),
                                        '../Metadata/merged_masking.bed')
    print(return_value)


def one_time_usage():
    return_value = check_regions_masked(Arm([Segment('Chr12', 64678139, 68251745)], 'test_arm'),
                                        '../Metadata/merged_masking.bed')
    print(return_value)


if __name__ == "__main__":
    test()
