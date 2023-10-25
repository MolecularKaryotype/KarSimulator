from Start_Genome import generate_genome_from_KT
from Structures import Segment


def kt_haploid_to_diploid(input_file, output_dir):
    """
    insert a 23X WT haploid to a modified KT (only modifying the KT section of the file)
    :param input_file:
    :param output_dir:
    :return:
    """
    genome = generate_genome_from_KT(input_file)
    if 'ChrX' in genome.full_KT:
        female_haploid = True
    else:
        female_haploid = False

    chr_of_interest = ['Chr' + str(i) for i in range(1, 23)]
    chr_of_interest.append('ChrX')
    if not female_haploid:
        # segment_dict_values = segment_dict.values()
        # max_dict_index = -1
        # for dict_index in segment_dict_values:
        #     try:
        #         value = int(dict_index)
        #         if value > max_dict_index:
        #             max_dict_index = value
        #     except ValueError:
        #         # Ignore non-integer elements
        #         pass
        # segment_dict[Segment('ChrX', 10000, 58605579)] = str(max_dict_index + 1)
        # segment_dict[Segment('ChrX', 58605580, 62412542)] = 'CENX'
        # segment_dict[Segment('ChrX', 62412543, 156030894)] = str(max_dict_index + 2)
        genome.motherboard.segments.append(Segment('ChrX', 10000, 58605579))
        genome.centromere_segments.append(Segment('ChrX', 58605580, 62412542))
        genome.motherboard.segments.append(Segment('ChrX', 62412543, 156030894))
        genome.motherboard.segments = sorted(genome.motherboard.segments)
    segment_dict = genome.segment_indexing()
    sorted_segments = sorted(segment_dict)

    return_str = 'chromosome\tKT\ttelo1_len\ttelo2_len\n'

    for current_chr in chr_of_interest:
        if current_chr == 'ChrX' and not female_haploid:
            haploid_chr_letter = 'a'
            t1_len = 10000
            t2_len = 10000
        else:
            haploid_chr_count = len(genome.full_KT[current_chr])
            haploid_chr_letter = str(chr(97 + haploid_chr_count))
            t1_len = -1
            t2_len = -1
            for genome_chr in genome.full_KT[current_chr]:
                t1_len = genome_chr.t1_len
                t2_len = genome_chr.t2_len
                if genome_chr.deleted:
                    return_str += '{}\tdeleted\t0\t0\n'.format(genome_chr.name)
                    continue
                tostring_segment_list = []
                for segment_itr in genome_chr:
                    if segment_itr.direction():
                        tostring_segment_list.append(segment_dict[segment_itr] + '+')
                    else:
                        new_segment_itr = segment_itr.duplicate()
                        new_segment_itr.invert()
                        tostring_segment_list.append(segment_dict[new_segment_itr] + '-')

                return_str += '{}\t{}\t{}\t{}\n'.format(genome_chr.name, ','.join(tostring_segment_list),
                                                        str(genome_chr.t1_len), str(genome_chr.t2_len))

        # insert WT haploid
        current_segments = []
        for segment_itr in sorted_segments:
            if segment_itr.chr_name == current_chr:
                current_segments.append(segment_itr)
        current_segment_indices = []
        for segment_itr in current_segments:
            current_segment_indices.append(segment_dict[segment_itr] + "+")

        return_str += '{}\t{}\t{}\t{}\n'.format(current_chr + haploid_chr_letter, ','.join(current_segment_indices),
                                                str(t1_len), str(t2_len))

    # add ChrY, which never gets the WT half
    if not female_haploid:
        for genome_chr in genome.full_KT['ChrY']:
            if genome_chr.deleted:
                return_str += '{}\tdeleted\t0\t0\n'.format(genome_chr.name)
                continue
            tostring_segment_list = []
            for segment_itr in genome_chr:
                if segment_itr.direction():
                    tostring_segment_list.append(segment_dict[segment_itr] + '+')
                else:
                    new_segment_itr = segment_itr.duplicate()
                    new_segment_itr.invert()
                    tostring_segment_list.append(segment_dict[new_segment_itr] + '-')

            return_str += '{}\t{}\t{}\t{}\n'.format(genome_chr.name, ','.join(tostring_segment_list),
                                                    str(genome_chr.t1_len), str(genome_chr.t2_len))

    file_raw_name = input_file.split('/')[-1]
    with open(output_dir + '/' + file_raw_name, 'w') as fp_write:
        fp_write.write(genome.motherboard_tostring())
        fp_write.write('---\n')
        fp_write.write(return_str)
        fp_write.write('---\n')
        fp_write.write(genome.history_tostring())


def cmd():
    import argparse
    parser = argparse.ArgumentParser(description="haploid KT to diploid KT")
    parser.add_argument("input_file")
    parser.add_argument("output_dir")
    args = parser.parse_args()

    kt_haploid_to_diploid(args.input_file, args.output_dir)


if __name__ == "__main__":
    cmd()
