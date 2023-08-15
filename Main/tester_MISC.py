from IO.read_FASTA import read_FASTA

seq_dict = read_FASTA('../Genomes/GCF_000001405.26_GRCh38_genomic.fasta',
                      ['NC_000001.11 Homo sapiens chromosome 1, GRCh38 Primary Assembly'])
print(seq_dict['NC_000001.11 Homo sapiens chromosome 1, GRCh38 Primary Assembly'][10000:10100])
