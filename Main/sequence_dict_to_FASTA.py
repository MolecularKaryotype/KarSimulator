def sequence_dict_to_FASTA(sequence_dict, output_file):
    with open(output_file, 'w') as fp_write:
        def custom_sort_chr(key):
            chr_part = key[3:-1]  # Extract the part after "Chr"
            if chr_part.isdigit():
                return int(chr_part)
            elif chr_part == "X":
                return int(23)  # Put ChrX at the end
            elif chr_part == "Y":
                return int(24)  # Put ChrY after ChrX
            return key
        sorted_keys = sorted(sequence_dict.keys(), key=custom_sort_chr)

        for header in sorted_keys:
            sequence = sequence_dict[header]
            fp_write.writelines(">{}\n".format(header))
            for i in range(0, len(sequence), 80):
                chunk = sequence[i:i+80]
                fp_write.writelines(chunk + "\n")
