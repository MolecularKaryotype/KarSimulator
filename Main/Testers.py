from Structures import *
from Start_Genome import generate_genome_from_KT
import IO


def test_KT_FASTA_Correspondence(KT_file, FASTA_file, genome_file):
    genome_KT = generate_genome_from_KT(KT_file)
    raw_genome = IO.read_FASTA(genome_file)
    pass
