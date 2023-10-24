from Structures import Segment, Arm
from read_masking_regions import read_masking_regions


def check_regions_masked(input_arm: Arm, masking_file):
    masking_arm = read_masking_regions(masking_file)
    if input_arm.arm_intersection(masking_arm):
        return input_arm.report_arm_intersection(masking_arm)
    else:
        return "Input Segments: " + str(input_arm) + "\nno intersection\n"


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


def cmd():
    import argparse
    parser = argparse.ArgumentParser(description="forbidden region intersection test")
    parser.add_argument("--segment", type=str, dest='input_segment', help="in format Chr1,10000,20000")
    args = parser.parse_args()

    masking_file = '../Metadata/merged_masking_unique.bed'
    input_parsed = args.input_segment.split(',')
    input_segment = Segment(input_parsed[0], int(input_parsed[1]), int(input_parsed[2]))
    print(check_regions_masked(Arm([input_segment], "check_masked"), masking_file))


def batch_check_masked():
    masking_file = '../Metadata/merged_masking_unique.bed'

    # 22q11.2 distal deletion syndrome
    input_segment = Segment('Chr22', 21562828, 23380258)
    print(check_regions_masked(Arm([input_segment], "check_masked"), masking_file))

    # 22q11 duplication syndrome
    input_segment = Segment('Chr22', 19022279, 21098156)
    print(check_regions_masked(Arm([input_segment], "check_masked"), masking_file))

    # 15q26 overgrowth syndrome
    input_segment = Segment('Chr15', 98814741, 101981189)
    print(check_regions_masked(Arm([input_segment], "check_masked"), masking_file))

    # 2p15-16.1 microdeletion syndrome
    input_segment = Segment('Chr2', 59058561, 61592680)
    print(check_regions_masked(Arm([input_segment], "check_masked"), masking_file))

    # Cri du Chat Syndrome (5p deletion)
    input_segment = Segment('Chr5', 10001, 12533192)
    print(check_regions_masked(Arm([input_segment], "check_masked"), masking_file))

    # WAGR 11p13 deletion syndrome
    input_segment = Segment('Chr11', 31784791, 32435541)
    print(check_regions_masked(Arm([input_segment], "check_masked"), masking_file))


def check_sunnyside():
    masking_file = '../Metadata/merged_masking_unique.bed'

    # 204
    input_segment = Segment('Chr15', 19024883, 22276055)
    print(check_regions_masked(Arm([input_segment], "check_masked"), masking_file))

    # 219
    input_segment = Segment('Chr11', 48350254, 48949944)
    print(check_regions_masked(Arm([input_segment], "check_masked"), masking_file))

    # 1327
    input_segment = Segment('ChrY', 2465162, 16252722)
    print(check_regions_masked(Arm([input_segment], "check_masked"), masking_file))
    input_segment = Segment('ChrY', 16264213, 56951882)
    print(check_regions_masked(Arm([input_segment], "check_masked"), masking_file))

    # 1338
    input_segment = Segment('Chr10', 34883020, 35028718)
    print(check_regions_masked(Arm([input_segment], "check_masked"), masking_file))


def check_keyhole():
    masking_file = '../Metadata/merged_masking_unique.bed'

    # 1878
    input_segment = Segment('Chr9', 28188952, 28352275)
    print(check_regions_masked(Arm([input_segment], "check_masked"), masking_file))

    # 2436
    input_segment = Segment('Chr9', 28171542, 28364123)
    print(check_regions_masked(Arm([input_segment], "check_masked"), masking_file))


if __name__ == "__main__":
    check_keyhole()
