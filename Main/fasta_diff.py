from read_FASTA import read_FASTA


def fasta_diff(fasta1, fasta2, header_conversion: {str: str}):
    """
    For the header pairs mentioned in header_conversion, check if two fasta files contain the same sequences
    :param fasta1: file 1
    :param fasta2: file 2
    :param header_conversion: a dictionary that converts fasta1's sequence header into fasta2's
    :return: None; output to stdout for each pair of headers in header_conversion
    """
    fasta1_dict = read_FASTA(fasta1, list(header_conversion.keys()))
    # for key in fasta1_dict:
    #     print(key, len(fasta1_dict[key]))

    fasta2_dict = read_FASTA(fasta2, list(header_conversion.values()))
    # for key in fasta2_dict:
    #     print(key, len(fasta2_dict[key]))

    for key, value in header_conversion.items():
        print("{} compared to {}: {}".format(key, value, fasta1_dict[key] == fasta2_dict[value]))


def test():
    # Create a dictionary to map chromosome names to modified names
    chromosome_mapping = {}

    # Define the range of chromosome names you want to modify
    chromosomes_to_modify = list(range(1, 23)) + ['X']

    # Iterate through the range and add entries to the dictionary
    for chromosome in chromosomes_to_modify:
        original_name = f'Chr{chromosome}'
        modified_name = f'Chr{chromosome}a'
        chromosome_mapping[original_name] = modified_name

    fasta_diff("../Genomes/hg38.fasta", "/Users/zhaoyangjia/Downloads/FASTAs/23X_ref.fasta", chromosome_mapping)

if __name__ == "__main__":
    test()
