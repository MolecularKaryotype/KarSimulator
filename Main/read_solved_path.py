from Structures import *


def read_solved_path(file):
    segment_dict = {}
    with open(file) as fp_read:
        fp_read.readline()

        for line in fp_read:
            line = line.replace("\n", "").split('\t')

            # documenting segments
            if line[0] == "Segment":
                chr_name = str(line[2])
                if chr_name == "23":
                    chr_name = "X"
                elif chr_name == "24":
                    chr_name = "Y"
                chr_name = "Chr" + chr_name
                start = int(line[3].split(".")[0])
                end = int(line[4].split(".")[0])
                segment_dict[int(line[1])] = Segment(chr_name, start, end, "solved_path")
            elif line[0].startswith("Path"):
                line = line[0].split(" = ")[1]
                line = line[:-1].split(" ")
                print(line)


def test():
    read_solved_path("/media/zhaoyang-new/workspace/KarSim/OMKar_outputs/simulation_final/NF1_microdeletion_r1.1/NF1_microdeletion_r1.1.txt")


if __name__ == "__main__":
    test()
