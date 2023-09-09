def read_FASTA(genome_path: str, chr_of_interest: [str]):
    sequence_dict = {}
    with open(genome_path) as fp_read:
        recording = False
        for line in fp_read:
            if line[0] == '>':
                # header line
                if recording:
                    sequence_dict[header] = ''.join(segment_sequence)
                    recording = False
                if line[1:].replace('\n', '') in chr_of_interest or chr_of_interest[0].lower() == 'all':
                    recording = True
                    segment_sequence = []
                    header = line[1:].replace('\n', '')
            else:
                # sequence line
                if recording:
                    segment_sequence.append(line.replace('\n', ''))
    # record the last segment
    if recording:
        sequence_dict[header] = ''.join(segment_sequence)

    return sequence_dict


# d1 = read_FASTA('../Preparation/Chr1.fasta', ['Chr1'])
# d2 = read_FASTA('../Genomes/GCF_000001405.26_GRCh38_genomic.fasta',
#                 ['NC_000001.11 Homo sapiens chromosome 1, GRCh38 Primary Assembly'])
# print(d1['Chr1'] == d2['NC_000001.11 Homo sapiens chromosome 1, GRCh38 Primary Assembly'])

def test():
    all_chrom = ['Chr1', 'Chr2', 'Chr3', 'Chr4', 'Chr5', 'Chr6', 'Chr7', 'Chr8', 'Chr9', 'Chr10', 'Chr11',
                'Chr12', 'Chr13', 'Chr14', 'Chr15', 'Chr16', 'Chr17', 'Chr18', 'Chr19', 'Chr20', 'Chr21',
                'Chr22', 'ChrX', 'ChrY']
    seq_dict = read_FASTA("../Genomes/hg38.fasta", all_chrom)
    for chrom in all_chrom:
        # seq_dict['ChrX'][10000:2406768]
        for index, char in enumerate(seq_dict[chrom]):
            if char.lower() not in ['a', 'c', 'g', 't', 'n', 'm', 'r', 'k', 'w', 'y', 's', 'b', 'v', 'h', 'd']:
                print(chrom, index, char)


if __name__ == "__main__":
    test()
