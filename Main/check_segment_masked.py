from Structures import Segment, Arm
from read_masking_regions import read_masking_regions


def check_regions_masked(input_arm: Arm, masking_file):
    masking_arm = read_masking_regions(masking_file)
    if input_arm.arm_intersection(masking_arm):
        return input_arm.report_arm_intersection(masking_arm)
    else:
        return "no intersection"


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
    parser.add_argument("--segment", type=str, dest='input_segment', help="in format (Chr1,10000,20000)")
    args = parser.parse_args()

    masking_file = '../Metadata/merged_masking_unique.bed'
    total_input = args.input_segment.split('(')
    for segment_index in range(1, len(total_input)):
        input_parsed = total_input[segment_index]
        input_parsed = input_parsed.replace('(', '').replace(')', '').replace(',', '').split('-')
        input_segment = Segment(input_parsed[0], int(input_parsed[1]), int(input_parsed[2]))
        print(check_regions_masked(Arm([input_segment], "check_masked"), masking_file))


def batch_check_masked_1020():
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


def batch_check_masked_1025():
    masking_file = '../Metadata/merged_masking_unique.bed'

    print('CMT1A')
    input_segment = Segment('Chr17', 14194598, 15567589)
    print(check_regions_masked(Arm([input_segment], "check_masked"), masking_file))

    print('NF1')
    input_segment = Segment('Chr17', 30780079, 31936302)
    print(check_regions_masked(Arm([input_segment], "check_masked"), masking_file))

    print('Potocki-Lupski')
    input_segment = Segment('Chr17', 16869758, 20318836)
    print(check_regions_masked(Arm([input_segment], "check_masked"), masking_file))

    print('Potocki-Shaffer')
    input_segment = Segment('Chr11', 43973250, 46030899)
    print(check_regions_masked(Arm([input_segment], "check_masked"), masking_file))

    print('STS')
    input_segment = Segment('ChrX', 6537771, 8165154)
    print(check_regions_masked(Arm([input_segment], "check_masked"), masking_file))

    print('RCAD')
    input_segment = Segment('Chr17', 36459259, 37856298)
    print(check_regions_masked(Arm([input_segment], "check_masked"), masking_file))

    print('Cri_du_chat')
    input_segment = Segment('Chr5', 10001, 12533192)
    print(check_regions_masked(Arm([input_segment], "check_masked"), masking_file))

    print('Angelman')
    input_segment = Segment('Chr15', 22677345, 28193120)
    print(check_regions_masked(Arm([input_segment], "check_masked"), masking_file))


def batch_check_masked_ddd():
    masking_file = '../Metadata/merged_masking_unique.bed'

    print('\n26674')
    input_segment = Segment('Chr8', 7167369, 8173892)
    print(input_segment)
    print(check_regions_masked(Arm([input_segment], "check_masked"), masking_file))

    print('\n43992')
    input_segment = Segment('Chr15', 19939373, 21061292)
    print(input_segment)
    print(check_regions_masked(Arm([input_segment], "check_masked"), masking_file))

    print('\n47063')
    input_segment = Segment('Chr16', 32117844, 32893189)
    print(input_segment)
    print(check_regions_masked(Arm([input_segment], "check_masked"), masking_file))

    print('\n42082')
    input_segment = Segment('Chr14', 19274726, 19942291)
    print(input_segment)
    print(check_regions_masked(Arm([input_segment], "check_masked"), masking_file))

    print('\n43992')
    input_segment = Segment('Chr15', 19939373, 21061292)
    print(input_segment)
    print(check_regions_masked(Arm([input_segment], "check_masked"), masking_file))

    print('\n26674')
    input_segment = Segment('Chr8', 7167369, 8173892)
    print(input_segment)
    print(check_regions_masked(Arm([input_segment], "check_masked"), masking_file))

    print('\n56143')
    input_segment = Segment('Chr22', 18170680, 18809686)
    print(input_segment)
    print(check_regions_masked(Arm([input_segment], "check_masked"), masking_file))

    print('\n47063')
    input_segment = Segment('Chr16', 32117844, 32893189)
    print(input_segment)
    print(check_regions_masked(Arm([input_segment], "check_masked"), masking_file))


def batch_check_masked_1026():
    masking_file = '../Metadata/merged_masking_unique.bed'

    print('\n8q21.11 Microdeletion Syndrome')
    input_segment = Segment('Chr8', 76314229, 76854003)
    print(input_segment.thousand_delimited())
    print(check_regions_masked(Arm([input_segment], "check_masked"), masking_file))

    print('\nEarly-onset Alzheimer disease with cerebral amyloid angiopathy')
    input_segment = Segment('Chr21', 25880549, 26171128)
    print(input_segment.thousand_delimited())
    print(check_regions_masked(Arm([input_segment], "check_masked"), masking_file))

    print('\nXp11.22-p11.23 Microduplication')
    input_segment = Segment('ChrX', 48476161, 52374518)
    print(input_segment.thousand_delimited())
    print(check_regions_masked(Arm([input_segment], "check_masked"), masking_file))

    print('\n1q21.1 recurrent microdeletion (susceptibility locus for neurodevelopmental disorders)')
    input_segment = Segment('Chr1', 147061832, 148411223)
    print(input_segment.thousand_delimited())
    print(check_regions_masked(Arm([input_segment], "check_masked"), masking_file))

    print('\n12q14 microdeletion syndrome')
    input_segment = Segment('Chr12', 64678139, 68251745)
    print(input_segment.thousand_delimited())
    print(check_regions_masked(Arm([input_segment], "check_masked"), masking_file))

    print('\n1q21.1 recurrent microdeletion (susceptibility locus for neurodevelopmental disorders)')
    input_segment = Segment('Chr1', 147061832, 148411223)
    print(input_segment.thousand_delimited())
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
    batch_check_masked_1026()
