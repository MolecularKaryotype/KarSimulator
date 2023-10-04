import IO


def Get_Masking_Regions(genome_path, output_path):
    # TODO: add filter against telomere and centromere
    sequence_dict = IO.read_FASTA(genome_path, ['all'])
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
                    # output
                    output_dict[itr_header].append(tuple([substring_start_index, current_index - 1]))
                    N_string = False
            current_index += 1

        if N_string:
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
    Get_Masking_Regions('','')

if __name__ == "__main__":
    test()



