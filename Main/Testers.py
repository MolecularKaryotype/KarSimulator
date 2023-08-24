from Structures import *
from Start_Genome import generate_genome_from_KT
import IO
import difflib


def test_KT_FASTA_Correspondence(KT_file, FASTA_file, genome_file):
    genome_KT = generate_genome_from_KT(KT_file)
    genome_FASTA = IO.read_FASTA(FASTA_file, ['all'])
    raw_genome = IO.read_FASTA(genome_file, ['all'])

    for chromosome in genome_KT:
        if chromosome.deleted:
            continue
        chromosome_sequence = ['N' * chromosome.t1_len]
        for segment in chromosome:
            if segment.direction():
                segment_sequence = raw_genome[segment.chr_name][segment.start: segment.end + 1]
            else:
                segment_sequence = raw_genome[segment.chr_name][segment.end: segment.start + 1][::-1]
            chromosome_sequence.append(segment_sequence)
        chromosome_sequence.append('N' * chromosome.t2_len)
        chromosome_sequence = ''.join(chromosome_sequence)

        print(chromosome.name)
        FASTA_correspondence_chromosome_sequence = genome_FASTA[chromosome.name]
        print(chromosome_sequence == FASTA_correspondence_chromosome_sequence)
        # for i, s in enumerate(difflib.ndiff(chromosome_sequence, FASTA_correspondence_chromosome_sequence)):
        #     if s[0] == ' ':
        #         continue
        #     elif s[0] == '-':
        #         print(u'Delete "{}" from position {}'.format(s[-1], i))
        #     elif s[0] == '+':
        #         print(u'Add "{}" to position {}'.format(s[-1], i))

