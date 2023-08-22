# KarSimulator

### Version v0.1

## Table of contents:
1. [Introduction](#Introduction)
2. [Workflow](#Workflow)
3. [Supported SVs](#Supported SVs)
4. [Installation](#Installation)
5. [Prerequisites](#Prerequisites)
6. [Usage](#Usage)

## Introduction
This tool randomly generates structural variations from a given genome, outputting the karyotype of the 
modified genome and a comprehensive event history. 

## Workflow
- Produce an unedited genome from a FASTA file, or load one of the pre-edited hg38 disorder 
as the starting genome.
- Edit the starting genome using a set of random parameters. 
This process yields `m` copies of diverse random karyotypes, 
all under the same set of random parameters. Each karyotype contains `n` structural variations.
- Choose one of the output karyotypes as the starting point for the next round of random generation. 
Apply a second set of random parameters. A comprehensive history that follow through both edits will be outputted.
- Repeat the process of generating layers of events. 
Then, select the karyotype of interest and output the corresponding FASTA.

## Supported SVs
![Supported SVs](/pics/Supported_SV.png)

## Installation
`git clone https://github.com/Zhaoyang-Jia/KarSimulator.git`

## Prerequisites
1. Python 3.9+
2. Genome file deposited in `Genomes/`
3. Genome Index file deposited in `Genomes/`
4. Dgenies (Optional for alignment plot)

## Usage
### Prepare a starting genome
- Load a prepared starting genome, this includes standard male and female and some common disorders:
   `python KarSimulator/src/KarSimulator.py selectGenome`
- Start with a custom unedited genome: `python KarSimulator/src/KarSimulator.py rawGenome`

#### Arguments:
| Argument  | Type  | Description                                                                                               |
|-----------|-------|-----------------------------------------------------------------------------------------------------------|
| `--copy`  | STR   | Copy number for the autosomes                                                                             |
| `--auto`  | [STR] | `[all]` for all 22 autosomes, `[Chr1, Chr2, etc.]` for custom selection                                   |
| `--sex`   | [STR] | `[male]`, `[female]` for XY and XX,  `[ChrX, ChrY, etc.]` for custom selection                            |
| `--index` | FILE  | (Optional) Genome Index File, only use when using non-hg38 genome <br/>(default: pre-compiled hg38 index) |

### Work on an input genome
- Full random mode: `python KarSimulator/src/KarSimulator.py random`
#### Arguments:
| Argument   | Type | Description                                                        |
|------------|------|--------------------------------------------------------------------|
| `--json`   | FILE | JSON file containing Random Mode parameters                        |
| `--genome` | FILE | Karyotype file containing the input karyotype                      |
| `-o`       | PATH | (Optional) output directory (default: same directory as JSON file) |
**NOTE1:** template JSON file can be found in `Sample_JSON/`
- Manual mode: not recommended in later layers as the indices are relative to the current genome, not the 
standard hg38 genome. Manual calculation needs to be done to account for off-shifts.
(e.g. if a deletion is to be made on Chr1:100,000-200,000 on a genome that had a prior deletion
on Chr1:50,000-100,000, the manual parameter needs to address deletion range Chr1:50,000-150,000):
`python KarSimulator/src/KarSimulator.py manual`
#### Arguments:
| Argument        | Type  | Description                                                                |
|-----------------|-------|----------------------------------------------------------------------------|
| `--instuctions` | FILE  | Python script file containing Manual mode scripts, to be executed in-order |
| `--genome`      | FILE  | Karyotype file containing the input karyotype                              |

#### Instructions:
Format in development

### Output FASTA
`python KarSimulator/src/KarSimulator.py fasta`
#### Arguments:
| Argument   | Type  | Description                                                             |
|------------|-------|-------------------------------------------------------------------------|
| `--genome` | FILE  | Karyotype file containing the input karyotype                           |
| `-o`       | FILE  | (Optional) output directory (default: same directory as Karyotype file) |
