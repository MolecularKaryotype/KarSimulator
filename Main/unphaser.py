from Start_Genome import generate_genome_from_KT
from Start_Genome import generate_raw_genome


def unphaser(kt_file):
    genome = generate_genome_from_KT(kt_file)

    ### extract initialization and generate new KT
    initialization_blocks = genome.initialization_string.split('\n')
    initialization_autosomes = initialization_blocks[1].replace('autosomes: ', '').replace(' ', '').replace('\t', '')
    initialization_copy_number = int(initialization_blocks[2].replace('\tautosomal copy number: ', ''))
    initialization_sex_chromosomes = initialization_blocks[3].replace('\tsex chromosomes: ', '').replace(' ', '')

    initialization_autosomes = initialization_autosomes.replace('[', '').replace(']', '').replace("\'", '').split(',')
    initialization_sex_chromosomes = initialization_sex_chromosomes.replace('[', '').replace(']', '').replace("\'",
                                                                                                              '').split(
        ',')
    new_genome = generate_raw_genome(initialization_copy_number, initialization_autosomes,
                                     initialization_sex_chromosomes, '../Genomes/hg38_index.txt')

    ### get dependent chromosome clusters that require manual unphasing
    history_segment_chr = {}
    for history_itr in genome.history:
        for history_segment_itr in history_itr[1].segments:
            if history_segment_itr in history_segment_chr:
                history_segment_chr[history_segment_itr].append(history_itr[2].name)
            else:
                history_segment_chr[history_segment_itr] = [history_itr[2].name]
            history_segment_chr[history_segment_itr].append(history_itr[3].name)

    segment_to_skip = []
    for segment_itr in history_segment_chr:
        current_chr_list = history_segment_chr[segment_itr]
        end_iteration = False
        for chr1_index in range(len(current_chr_list)):
            for chr2_index in range(chr1_index + 1, len(current_chr_list)):
                chr1_number = current_chr_list[chr1_index][:-1]
                chr2_number = current_chr_list[chr2_index][:-1]
                if chr1_number == chr2_number:
                    chr1_homolog = current_chr_list[chr1_index][-1]
                    chr2_homolog = current_chr_list[chr2_index][-1]
                    if chr1_homolog != chr2_homolog:
                        segment_to_skip.append(segment_itr)
                        end_iteration = True
                        break

            if end_iteration:
                break

    # all chr involving the segment_to_skip are dependent and pending manual review
    chr_to_skip = set()
    for segment_itr in segment_to_skip:
        for chr_itr in history_segment_chr[segment_itr]:
            chr_to_skip.add(chr_itr[:-1])

    ### perform automatic unphasing on independent homolog segments
    # TODO: how to keep track of the relative location of the events??? translocation is unsolvable
    for history_itr in genome.history:
        event_name = history_itr[0]

        if event_name == 'deletion':
            pass


def test():
    unphaser("../test_folder/1q21-1_recurrent_microdeletion_v2_r1.kt.txt")


if __name__ == "__main__":
    test()
