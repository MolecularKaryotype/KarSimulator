def sequence_dict_to_FASTA(sequence_dict, output_file):
    with open(output_file, 'w') as fp_write:
        for header, sequence in sequence_dict.items():
            fp_write.writelines(">{}\n".format(header))
            for i in range(0, len(sequence), 80):
                chunk = sequence[i:i+80]
                fp_write.writelines(chunk + "\n")
