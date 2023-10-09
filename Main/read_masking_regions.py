from Main.Structures import Arm, Segment


def read_masking_regions(masking_file) -> Arm:
    segment_list = []
    with open(masking_file) as fp_read:
        fp_read.readline()
        for line in fp_read:
            line = line.replace('\n', '').split('\t')
            line[1] = line[1].split('.')[0]  # remove the .5
            line[2] = line[2].split('.')[0]
            new_segment = Segment('Chr' + str(line[0]), int(line[1]), int(line[2]))
            segment_list.append(new_segment)
    return Arm(segment_list, 'masking_regions')


def test():
    print(read_masking_regions('../Metadata/test_masking_regions.bed'))


if __name__ == "__main__":
    test()
