from Main.read_FASTA import read_FASTA


def Get_Masking_Regions(genome_path, output_path):
    # TODO: add filter against telomere and centromere
    sequence_dict = read_FASTA(genome_path, ['all'])
    output_dict = {key: [] for key in sequence_dict.keys()}

    for itr_header in sequence_dict:
        sequence = sequence_dict[itr_header]
        current_index = 0
        substring_start_index = -1

        N_string = False
        while current_index < len(sequence):
            if not N_string and sequence[current_index] == 'N':
                N_string = True
                substring_start_index = current_index
            elif N_string:
                if sequence[current_index] != 'N':
                    # output, ignore regions that are too small
                    if current_index - substring_start_index >= 50000:
                        output_dict[itr_header].append(tuple([substring_start_index, current_index - 1]))
                    N_string = False
            current_index += 1

        if N_string:
            if len(sequence) - substring_start_index >= 50000:
                output_dict[itr_header].append(tuple([substring_start_index, len(sequence) - 1]))

    with open(output_path, 'w') as fp_write:
        for key in output_dict:
            value = output_dict[key]
            fp_write.write(str(key) + ':\n')
            for pair in value:
                fp_write.write(str(pair) + '\n')


def test():
    Get_Masking_Regions('test_Get_Masking_Regions', 'test_output')


def generate_hg38():
    Get_Masking_Regions('/media/zhaoyang-new/workspace/KarSim/KarSimulator/Genomes/hg38.fasta',
                        '../Metadata/hg38_masking_regions_50000.txt')


if __name__ == "__main__":
    generate_hg38()



