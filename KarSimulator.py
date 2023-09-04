import argparse
import json
import random
import os
import shutil

from Main.Start_Genome import *


def rawGenome_mode(args):
    print("Running rawGenome mode with arguments:", args)
    autosome_list = args.autosomes.replace('[', '').replace(']', '').split(',')
    sex_chromosome_list = args.sex_chromosomes.replace('[', '').replace(']', '').split(',')
    genome = generate_raw_genome(args.copy_number, autosome_list, sex_chromosome_list, args.index_file)
    genome.output_KT(args.output_dir + args.name + '.kt.txt')


def random_mode(args):
    num_supported_SV = 12
    error_logs_path = "./error_logs/"
    print("Running random mode with arguments:", args)

    with open(args.json_file) as file:
        instruction = json.load(file)
    job_name = instruction['job_name']
    KT_input_file = instruction['template_KT']
    output_file_prefix = instruction['output_directory'] + instruction['output_file_name']
    event_settings = instruction['event_setting']
    number_of_events = instruction['number_of_events']
    number_of_iterations = instruction['number_of_iterations']

    # scale event weights
    sum_weight = 0.0
    for index in range(num_supported_SV):
        sum_weight += event_settings[index]['likelihood_weight']
    event_likelihoods = []
    for index in range(num_supported_SV):
        event_likelihoods.append(float(event_settings[index]['likelihood_weight']) / sum_weight)

    for event_iteration_index in range(number_of_iterations):
        # dump error logs
        for file in os.listdir(error_logs_path):
            file_path = os.path.join(error_logs_path, file)
            if os.path.isfile(file_path):
                os.remove(file_path)
        # copy KT input into error logs
        shutil.copy(KT_input_file, error_logs_path)
        random_parameter_error_logs = os.path.join(error_logs_path, 'random_parameters.txt')
        open(random_parameter_error_logs, 'w').close()

        def error_log_tostring(function_name, *function_args):
            params = ', '.join(map(str, function_args))
            return f"{function_name}: {params}"

        genome = generate_genome_from_KT(KT_input_file)
        full_output_file_path = output_file_prefix + "_r" + str(event_iteration_index + 1) + ".kt.txt"
        for event_index in range(number_of_events):
            # choose event
            current_event = random.choices(range(num_supported_SV), weights=event_likelihoods)[0]

            # choose chr
            chr_weights = []
            sum_length = 0.0
            for chromosome in genome:
                sum_length += len(chromosome)
            for chromosome in genome:
                chr_weights.append(float(len(chromosome)) / sum_length)
            current_chr1 = random.choices(genome.get_chromosome_list(), chr_weights)[0]
            current_chr2 = current_chr1
            if current_event != 6:
                while current_chr1 == current_chr2:
                    current_chr2 = random.choices(genome.get_chromosome_list(), chr_weights)[0]

            # choose arm
            arm1_weights = [float(current_chr1.p_arm_len()) / len(current_chr1),
                            float(current_chr1.q_arm_len()) / len(current_chr1)]
            arm2_weights = [float(current_chr2.p_arm_len()) / len(current_chr2),
                            float(current_chr2.q_arm_len()) / len(current_chr2)]
            current_arm1_value = random.choices(["p", "q"], arm1_weights)[0]
            current_arm2_value = random.choices(["p", "q"], arm2_weights)[0]
            if current_arm1_value == 'p':
                current_arm1 = current_chr1.p_arm
            else:
                current_arm1 = current_chr1.q_arm
            if current_arm2_value == 'p':
                current_arm2 = current_chr2.p_arm
            else:
                current_arm2 = current_chr2.q_arm

            # choose length
            current_event1_length = -1
            current_event2_length = -1
            if current_event not in [7, 8, 9, 10, 11]:
                current_event1_length = random.randint(event_settings[current_event]['min_size'],
                                                       event_settings[current_event]['max_size'])
                current_event1_length = min(current_event1_length, len(current_arm1) - 1)
                if current_event in [5, 6]:
                    current_event2_length = random.randint(event_settings[current_event]['min_size2'],
                                                           event_settings[current_event]['max_size2'])
                    current_event2_length = min(current_event2_length, len(current_arm2) - 1)

            # choose start location
            current_event_start_location1 = -1
            current_event_start_location2 = -1
            if current_event not in [7, 8, 9, 10, 11]:
                current_event_start_location1 = random.randint(0, len(current_arm1) - current_event1_length - 1)
                if current_event in [5, 6]:
                    current_event_start_location2 = random.randint(0, len(current_arm2) - current_event2_length - 1)

            def run_terminal_likelihood():
                terminal_likelihood = event_settings[current_event]['terminal_occurrence_likelihood']
                non_terminal_likelihood = 1 - terminal_likelihood
                terminal_status = \
                    random.choices([True, False], [terminal_likelihood, non_terminal_likelihood])[0]
                return terminal_status

            # perform event
            if current_event == 0:
                terminal_event = run_terminal_likelihood()
                if terminal_event:
                    if current_arm1 == current_chr1.p_arm:
                        current_event_start_location1 = 0
                    else:
                        current_event_start_location1 = len(current_arm1) - current_event1_length

                with open(random_parameter_error_logs, 'a') as fp_write:
                    command_append = error_log_tostring('deletion', current_arm1,
                                                        current_event_start_location1,
                                                        current_event_start_location1 + current_event1_length)
                    fp_write.write(command_append + '\n')
                genome.deletion(current_chr1, current_arm1,
                                current_event_start_location1,
                                current_event_start_location1 + current_event1_length)
            elif current_event == 1:
                terminal_event = run_terminal_likelihood()
                if terminal_event:
                    if current_arm1 == current_chr1.p_arm:
                        current_event_start_location1 = 0
                    else:
                        current_event_start_location1 = len(current_arm1) - current_event1_length

                with open(random_parameter_error_logs, 'a') as fp_write:
                    command_append = error_log_tostring('inversion', current_arm1,
                                                        current_event_start_location1,
                                                        current_event_start_location1 + current_event1_length)
                    fp_write.write(command_append + '\n')
                genome.inversion(current_chr1, current_arm1,
                                 current_event_start_location1,
                                 current_event_start_location1 + current_event1_length)
            elif current_event == 2:
                terminal_event = run_terminal_likelihood()
                if terminal_event:
                    if current_arm1 == current_chr1.p_arm:
                        current_event_start_location1 = 0
                    else:
                        current_event_start_location1 = len(current_arm1) - current_event1_length

                with open(random_parameter_error_logs, 'a') as fp_write:
                    command_append = error_log_tostring('tandem duplication', current_arm1,
                                                        current_event_start_location1,
                                                        current_event_start_location1 + current_event1_length)
                    fp_write.write(command_append + '\n')
                genome.duplication(current_chr1, current_arm1,
                                   current_event_start_location1,
                                   current_event_start_location1 + current_event1_length)
            elif current_event == 3:
                terminal_event = run_terminal_likelihood()
                if terminal_event:
                    if current_arm1 == current_chr1.p_arm:
                        current_event_start_location1 = 0
                    else:
                        current_event_start_location1 = len(current_arm1) - current_event1_length

                with open(random_parameter_error_logs, 'a') as fp_write:
                    command_append = error_log_tostring('segmental duplication', current_arm1,
                                                        current_event_start_location1,
                                                        current_event_start_location1 + current_event1_length)
                    fp_write.write(command_append + '\n')
                genome.duplication(current_chr1, current_arm1,
                                   current_event_start_location1,
                                   current_event_start_location1 + current_event1_length)
            elif current_event == 4:
                terminal_event = run_terminal_likelihood()
                if terminal_event:
                    if current_arm1 == current_chr1.p_arm:
                        current_event_start_location1 = 0
                    else:
                        current_event_start_location1 = len(current_arm1) - current_event1_length
                left_dup_inv_likelihood = event_settings[current_event]['left_dupinv_to_right_dupinv_likelihood']
                right_dup_inv_likelihood = 1 - left_dup_inv_likelihood
                event_direction = \
                    random.choices(['left', 'right'], [left_dup_inv_likelihood, right_dup_inv_likelihood])[0]

                if event_direction == 'left':
                    with open(random_parameter_error_logs, 'a') as fp_write:
                        command_append = error_log_tostring('left_duplication_inversion', current_arm1,
                                                            current_event_start_location1,
                                                            current_event_start_location1 + current_event1_length)
                        fp_write.write(command_append + '\n')
                    genome.left_duplication_inversion(current_chr1, current_arm1,
                                                      current_event_start_location1,
                                                      current_event_start_location1 + current_event1_length)
                else:
                    with open(random_parameter_error_logs, 'a') as fp_write:
                        command_append = error_log_tostring('right_duplication_inversion', current_arm1,
                                                            current_event_start_location1,
                                                            current_event_start_location1 + current_event1_length)
                        fp_write.write(command_append + '\n')
                    genome.right_duplication_inversion(current_chr1, current_arm1,
                                                       current_event_start_location1,
                                                       current_event_start_location1 + current_event1_length)
            elif current_event in [5, 6]:
                with open(random_parameter_error_logs, 'a') as fp_write:
                    command_append = error_log_tostring('reciprocal translocation', current_chr1, current_arm1,
                                                        current_event_start_location1,
                                                        current_event_start_location1 + current_event1_length,
                                                        current_chr2, current_arm2, current_event_start_location2,
                                                        current_event_start_location2 + current_event2_length)
                    fp_write.write(command_append + '\n')
                genome.translocation_reciprocal(current_chr1, current_arm1, current_event_start_location1,
                                                current_event_start_location1 + current_event1_length,
                                                current_chr2, current_arm2, current_event_start_location2,
                                                current_event_start_location2 + current_event2_length)
            elif current_event == 7:
                with open(random_parameter_error_logs, 'a') as fp_write:
                    command_append = error_log_tostring('arm_deletion', current_chr1, current_arm1)
                    fp_write.write(command_append + '\n')
                genome.arm_deletion(current_chr1, current_arm1)
            elif current_event in [8, 9]:
                with open(random_parameter_error_logs, 'a') as fp_write:
                    command_append = error_log_tostring('arm_tandem_duplication', current_chr1, current_arm1)
                    fp_write.write(command_append + '\n')
                genome.arm_tandem_duplication(current_chr1, current_arm1)
            elif current_event == 10:
                with open(random_parameter_error_logs, 'a') as fp_write:
                    command_append = error_log_tostring('chromosomal_deletion', current_chr1)
                    fp_write.write(command_append + '\n')
                genome.chromosomal_deletion(current_chr1)
            elif current_event == 11:
                with open(random_parameter_error_logs, 'a') as fp_write:
                    command_append = error_log_tostring('chromosomal_duplication', current_chr1)
                    fp_write.write(command_append + '\n')
                genome.chromosomal_duplication(current_chr1)

        genome.mark_history(job_name)
        genome.output_KT(full_output_file_path)


def manual_mode(args):
    print("Running manual mode with arguments:", args)
    print("currently under development")


def fasta_mode(args):
    if args.name is None:
        args.name = args.input_kar_file.split('/')[-1].split('.')[0]
    if args.output_dir is None:
        args.output_dir = '/'.join(args.input_kar_file.split('/')[:-1]) + '/'
    print("Running fasta mode with arguments:", args)
    genome = generate_genome_from_KT(args.input_kar_file)
    genome.output_FASTA(args.genome_file, args.output_dir + args.name + '.fasta')


def main():
    parser = argparse.ArgumentParser(description="KarSimulator command-line interface")
    subparsers = parser.add_subparsers(dest="mode", help="Choose a mode")

    # rawGenome mode
    rawGenome_parser = subparsers.add_parser("rawGenome", help="Run rawGenome mode: generate an unedited karyotype")
    rawGenome_parser.add_argument("--name", type=str, default='unnamed', dest="name",
                                  help="(default: unnamed) Name of the output KT file")
    rawGenome_parser.add_argument("--copy", type=int, default=2, dest="copy_number",
                                  help="(default: 2) Copy number for the autosomes")
    rawGenome_parser.add_argument("--auto", type=str, default='[all]', dest="autosomes",
                                  help="(default: all 22 autosomes) Autosomes selection, "
                                       "[all] for all 22 autosomes; [Chr1,Chr2,etc.] for custom selection")
    rawGenome_parser.add_argument("--sex", type=str, default='[female]', dest="sex_chromosomes",
                                  help="Sex chromosome selection, "
                                       "(default: female) [male],[female] for XY and XX, respectively; "
                                       "[ChrX,ChrY,etc.] for custom selection")
    rawGenome_parser.add_argument("--index", type=str, default='Genomes/hg38_index.txt', dest="index_file",
                                  help="(default: pre-compiled hg38 index) Genome Index File path")
    rawGenome_parser.add_argument("-o", type=str, default='./', dest="output_dir",
                                  help="(default: current directory) output directory")

    # random mode
    random_parser = subparsers.add_parser("random", help="Run random mode: introduces a set number of SVs on top of "
                                                         "current karyotype")
    random_parser.add_argument("--json", type=str, dest='json_file',
                               help="JSON file containing Random Mode parameters")

    # manual mode
    manual_parser = subparsers.add_parser("manual", help="Run manual mode: introduces a series of SVs on top of "
                                                         "current karyotype")
    manual_parser.add_argument("--json", type=str, dest='json_file',
                               help="JSON file containing Random Mode parameters")

    # fasta mode
    fasta_parser = subparsers.add_parser("fasta", help="Run fasta mode: output fasta from karyotype")
    fasta_parser.add_argument("--name", type=str, dest="name",
                              help="(default: the prefix of the input kar file) Name of the output FASTA file")
    fasta_parser.add_argument("--genome", type=str, dest='genome_file', default='./Genomes/hg38.fasta',
                              help="(default: gh38 in ./Genomes/) Genome FASTA file")
    fasta_parser.add_argument("--kar", type=str, dest='input_kar_file',
                              help="Karyotype file containing the input karyotype")
    fasta_parser.add_argument("-o", type=str, dest="output_dir",
                              help="(default: same directory as Karyotype file) output directory")

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
