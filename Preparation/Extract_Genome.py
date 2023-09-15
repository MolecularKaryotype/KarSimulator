import IO


def Extract_Genome(chr_of_interest: [str], genome_path: str, chr_name_file: str, output_file: str):
    chr_name_conversion = {}
    full_name_list = []
    if chr_name_file == "none":
        full_name_list = chr_of_interest
    else:
        with open(chr_name_file) as fp_read:
            for line in fp_read:
                line = line.replace("\n", "").split("\t")
                if line[0] in chr_of_interest:
                    chr_name_conversion[line[1]] = line[0]
                    full_name_list.append(line[1])

    sequence_dict = IO.read_FASTA(genome_path, full_name_list)

    new_sequence_dict = {}
    for header, sequence in sequence_dict.items():
        if chr_name_file == "none":
            new_header = header
        else:
            new_header = chr_name_conversion[header]
        new_sequence_dict[new_header] = sequence

    IO.sequence_dict_to_FASTA(new_sequence_dict, output_file)


# all_chr = ['Chr1', 'Chr2', 'Chr3', 'Chr4', 'Chr5', 'Chr6', 'Chr7', 'Chr8', 'Chr9', 'Chr10', 'Chr11',
#            'Chr12', 'Chr13', 'Chr14', 'Chr15', 'Chr16', 'Chr17', 'Chr18', 'Chr19', 'Chr20', 'Chr21',
#            'Chr22', 'ChrX', 'ChrY']
# Extract_Genome(all_chr, '../Genomes/GCF_000001405.26_GRCh38_genomic.fasta', '../Metadata/Chr_Names.txt',
#                'All24Chr.fasta')

# Extract_Genome(['Chr14', 'Chr21', 'ChrX', 'ChrY'],
#                '../Genomes/GCF_000001405.26_GRCh38_genomic.fasta', '../Metadata/Chr_Names.txt',
#                'Chr14Chr21ChrXY.fasta')


def one_time_usage():
    Extract_Genome(['ChrXa'],
                   '/media/zhaoyang-new/workspace/KarSim/0910_test/FASTA/23Xe10_r1.fasta', "none",
                   '/media/zhaoyang-new/workspace/KarSim/0915_test/23Xe10_r1_ChrX.fasta')


if __name__ == "__main__":
    one_time_usage()
