import argparse
import os
import json
import random

from Main.Start_Genome import *


def rawGenome_mode(args):
    print("Running rawGenome mode with arguments:", args)
    autosome_list = args.autosomes.replace('[', '').replace(']', '').split(',')
    sex_chromosome_list = args.sex_chromosomes.replace('[', '').replace(']', '').split(',')
    genome = generate_raw_genome(args.copy_number, autosome_list, sex_chromosome_list, args.index_file)
    genome.output_KT(args.output_dir)


def random_mode(args):
    print("Running random mode with arguments:", args)
    with open(args.json_file) as file:
        instruction = json.load(file)

    event_settings = instruction['event_setting']
    number_of_events = instruction['number_of_events']
    number_of_iterations = instruction['number_of_parallel_universe']
    # scale event weights
    sum_weight = 0.0
    for index in range(2):
        sum_weight += event_settings[index]['likelihood_weight']
    event_likelihoods = []
    for index in range(2):
        event_likelihoods.append(float(event_settings[index]['likelihood_weight']) / sum_weight)

    for event_iteration_index in range(number_of_iterations):
        genome = generate_genome_from_KT(args.input_kar_file)
        for event_index in range(number_of_events):
            # choose event
            current_event = random.choices(range(2), weights=event_likelihoods)[0]

            # choose chr
            current_chr1 = -1
            current_chr2 = -1
            chr_weights = []
            sum_length = 0.0
            for chromosome in genome:
                sum_length += len(chromosome)
            for chromosome in genome:
                chr_weights.append(float(len(chromosome)) / sum_length)
            current_chr1 = random.choices(genome.get_chromosome_list(), chr_weights)[0]

            # choose arm
            arm_weights = [float(current_chr1.p_arm_len()) / len(current_chr1),
                           float(current_chr1.q_arm_len()) / len(current_chr1)]
            current_arm1_value = random.choices(["p", "q"], arm_weights)[0]
            if current_arm1_value == 'p':
                current_arm1 = current_chr1.p_arm
            else:
                current_arm1 = current_chr1.q_arm

            # choose length
            current_event1_length = random.randint(event_settings[current_event]['min_size'],
                                                   event_settings[current_event]['max_size'])
            current_event1_length = min(current_event1_length, len(current_arm1) - 1)


def manual_mode(args):
    print("Running manual mode with arguments:", args)
    print("currently under development")


def fasta_mode(args):
    print("Running fasta mode with arguments:", args)
    genome = generate_genome_from_KT(args.input_kar_file)
    genome.output_FASTA(args.genome, args.output_dir)


def main():
    current_directory = os.getcwd()
    parent_directory = os.path.dirname(current_directory)
    os.chdir(parent_directory)

    parser = argparse.ArgumentParser(description="KarSimulator command-line interface")
    subparsers = parser.add_subparsers(dest="mode", help="Choose a mode")

    # rawGenome mode
    rawGenome_parser = subparsers.add_parser("rawGenome", help="Run rawGenome mode: generate an unedited karyotype")
    rawGenome_parser.add_argument("--copy", type=int, default=2, dest="copy_number",
                                  help="Copy number for the autosomes")
    rawGenome_parser.add_argument("--auto", type=str, default='[all]', dest="autosomes",
                                  help="Autosomes selection, "
                                       "[all] for all 22 autosomes; [Chr1,Chr2,etc.] for custom selection")
    rawGenome_parser.add_argument("--sex", type=str, default='[female]', dest="sex_chromosomes",
                                  help="Sex chromosome selection, "
                                       "[male],[female] for XY and XX, respectively; "
                                       "[ChrX,ChrY,etc.] for custom selection")
    rawGenome_parser.add_argument("--index", type=str, default='./Genome/hg38_index.txt', dest="index_file",
                                  help="Genome Index File path")
    rawGenome_parser.add_argument("-o", type=str, default='./', dest="output_dir",
                                  help="output directory")

    # random mode
    random_parser = subparsers.add_parser("random", help="Run random mode: introduces a set number of SVs on top of "
                                                         "current karyotype")
    random_parser.add_argument("--json", type=str, dest='json_file',
                               help="JSON file containing Random Mode parameters")
    random_parser.add_argument("--kar", type=str, dest='input_kar_file',
                               help="Karyotype file containing the input karyotype")
    random_parser.add_argument("-o", type=str, default='./', dest="output_dir",
                               help="output directory")

    # manual mode
    manual_parser = subparsers.add_parser("manual", help="Run manual mode: introduces a series of SVs on top of "
                                                         "current karyotype")
    manual_parser.add_argument("--json", type=str, dest='json_file',
                               help="JSON file containing Random Mode parameters")
    manual_parser.add_argument("--kar", type=str, dest='input_kar_file',
                               help="Karyotype file containing the input karyotype")
    manual_parser.add_argument("-o", type=str, default='./', dest="output_dir",
                               help="output directory")

    # fasta mode
    fasta_parser = subparsers.add_parser("fasta", help="Run fasta mode: output fasta from karyotype")
    fasta_parser.add_argument("--genome", type=str, dest='genome_file', default='./Genomes/hg38.fasta',
                              help="Genome FASTA file")
    fasta_parser.add_argument("--kar", type=str, dest='input_kar_file',
                              help="Karyotype file containing the input karyotype")
    fasta_parser.add_argument("-o", type=str, default='./', dest="output_dir",
                              help="output directory")

    # call corresponding mode
    args = parser.parse_args()
    if args.mode == "rawGenome":
        rawGenome_mode(args)
    elif args.mode == "random":
        random_mode(args)
    elif args.mode == "manual":
        manual_mode(args)
    elif args.mode == "fasta":
        fasta_mode(args)


if __name__ == "__main__":
    main()
