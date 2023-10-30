def translate_global_index_to_arm_index(chr_name, start_index, end_index, genome_index_file):
    p_start = -1
    q_start = -1

    with open(genome_index_file) as fp_read:
        for line in fp_read:
            line = line.replace('\n', '').split('\t')
            if line[0] == chr_name:
                # TODO: add warning for intersecting centromere and telomere regions
                p_start = int(line[2])
                q_start = int(line[4])

    if start_index >= q_start:
        print("{} q_arm\nstart: {}\nend: {}\n".format(chr_name, str(start_index - q_start), str(end_index - q_start)))
    else:
        print("{} p_arm\nstart: {}\nend: {}\n".format(chr_name, str(start_index - p_start), str(end_index - p_start)))


def batch_manual():
    indexing_file = "../Metadata/Full_Genome_Indices.txt"
    translate_global_index_to_arm_index('Chr22', 21562828, 23380258, indexing_file)

    translate_global_index_to_arm_index('Chr22', 19022279, 21098156, indexing_file)

    translate_global_index_to_arm_index('Chr15', 98814741, 101981189, indexing_file)

    translate_global_index_to_arm_index('Chr2', 59058561, 61592680, indexing_file)

    translate_global_index_to_arm_index('Chr5', 10001, 12533192, indexing_file)

    translate_global_index_to_arm_index('Chr11', 31784791, 32435541, indexing_file)


def batch_1025():
    indexing_file = "../Metadata/Full_Genome_Indices.txt"

    print('CMT1A')
    translate_global_index_to_arm_index('Chr17', 14194598, 15567589, indexing_file)

    print('NF1')
    translate_global_index_to_arm_index('Chr17', 30780079, 31936302, indexing_file)

    print('Potocki-Shaffer')
    translate_global_index_to_arm_index('Chr11', 43973250, 46030899, indexing_file)

    print('STS')
    translate_global_index_to_arm_index('ChrX', 6537771, 8165154, indexing_file)

    print('Cri_du_chat')
    translate_global_index_to_arm_index('Chr5', 10001, 12533192, indexing_file)

    print('Angelman')
    translate_global_index_to_arm_index('Chr15', 22677345, 28193120, indexing_file)


def batch_1026():
    indexing_file = "../Metadata/Full_Genome_Indices.txt"

    print('Xp11.22-p11.23 Microduplication')
    translate_global_index_to_arm_index('ChrX', 48476161, 52374518, indexing_file)

    print('Early-onset Alzheimer')
    translate_global_index_to_arm_index('Chr21', 25880549, 26171128, indexing_file)

    print('1q21.1 recurrent microduplication')
    translate_global_index_to_arm_index('Chr1', 147061832, 148411223, indexing_file)

    print('1q21.1 recurrent microdeletion')
    translate_global_index_to_arm_index('Chr1', 147061832, 148411223, indexing_file)

    print('8q21.11 Microdeletion')
    translate_global_index_to_arm_index('Chr8', 76314229, 76854003, indexing_file)

    print('12q14 microdeletion syndrome')
    translate_global_index_to_arm_index('Chr12', 64678139, 68251745, indexing_file)




if __name__ == "__main__":
    batch_1026()
