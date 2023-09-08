import IO


def fasta_diff(fasta1, fasta2, header_conversion: {str: str}):
    """
    For the header pairs mentioned in header_conversion, check if two fasta files contain the same sequences
    :param fasta1: file 1
    :param fasta2: file 2
    :param header_conversion: a dictionary that converts fasta1's sequence header into fasta2's
    :return: None; output to stdout for each pair of headers in header_conversion
    """
    fasta1_dict = read_FASTA(fasta1, list(header_conversion.keys()))
    fasta2_dict = read_FASTA(fasta2, list(header_conversion.values()))

    for key, value in header_conversion.items():
        print("{} compared to {}: {}".format(key, value, fasta1_dict[key] == fasta2_dict[value]))


def test():
    fasta_diff("fasta1.fasta", "fasta2.fasta", {'Chr1': "Chr1a", 'Chr2': 'Chr2a'})


if __name__ == "__main__":
    test()
