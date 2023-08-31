# KarSimulator

### Version v0.1

## Table of contents:
1. [Known Issues](#known-issues)
2. [Introduction](#introduction)
3. [Workflow](#workflow)
4. [Supported SVs](#supported-svs)
5. [Installation](#installation)
6. [Prerequisites](#prerequisites)
7. [Usage](#usage)

## Known Issues
If a bug is met while running random generation, please DO NOT proceed further. There is a constantly updating error log folder.
Please copy the content of the folder and share it with me, so I can replicate the issue.
Thank you!
1. Manual mode under development
2. segmental duplication and arm-segmental duplication disabled
3. non-reciprocal translocation under development
4. scope for the random event under development (to limit events to only happen on a certain
set of chromosomes)

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
2. Genome file deposited in `Genomes/` (`gunzip ./Genomes/hg38.fasta.gz` to unzip prepared hg38 genome)
3. Genome Index file deposited in `Genomes/` (hg38 indexing file is included)
4. Dgenies (Optional for alignment plot)

## Usage
### Prepare an unedited starting genome
- If genome of interest is included, skip this step. List of precompiled genomes is in the next step.
- Start with a custom unedited genome: `python ./Main/KarSimulator.py rawGenome`

#### Arguments:
| Argument  | Type  | Description                                                                                                |
|-----------|-------|------------------------------------------------------------------------------------------------------------|
| `--name`  | STR   | (default: unnamed) <br /> name of the output KT file                                                       |
| `--copy`  | STR   | (default: 2) <br /> Copy number for the autosomes                                                          |
| `--auto`  | [STR] | (default: all 22 autosomes) <br /> `[all]` for all 22 autosomes, `[Chr1, Chr2, etc.]` for custom selection |
| `--sex`   | [STR] | (default: female) <br /> `[male]`, `[female]` for XY and XX,  `[ChrX, ChrY, etc.]` for custom selection    |
| `--index` | FILE  | (default: pre-compiled hg38 index) <br /> Genome Index File                                                |
| `-o`      | PATH  | (default: current directory) <br /> output directory                                                       |

### Work on an input genome
- Full random mode: `python ./Main/KarSimulator.py random`
#### Arguments:
| Argument | Type | Description                                   |
|----------|------|-----------------------------------------------|
| `--json` | FILE | JSON file containing Random Mode parameters   |

#### Pre-compiled Karyotypes:
| Genome Path                                   | Description                                                                       | Structural Variations (SVs)                                                                |
|-----------------------------------------------|-----------------------------------------------------------------------------------|--------------------------------------------------------------------------------------------|
| ./Precompiled_Kar/Female.kt.txt               | Typical female genotype                                                           | N/a                                                                                        |
| ./Precompiled_Kar/Male.kt.txt                 | Typical male genotype                                                             | N/a                                                                                        |
| ./Precompiled_Kar/Down_Trisomy21.kt.txt       | Extra Chr21                                                                       | chromosomal_duplication(Chr21a)                                                            |
| ./Precompiled_Kar/Down_Extra21q.kt.txt        | Robertsonian Translocation of Chr21q onto Chr14q, resulting in 3 copies of Chr21q | duplication(Chr21a-q); deletion(Chr14a-q); nonreciprocal_translocation(Chr21a-q, Chr14a-q) |
| ./Precompiled_Kar/Prader_Willi.kt.txt         | Loss of Chr15q11-13                                                               | deletion(Chr15a:2050001-33400000)                                                          |
| ./Precompiled_Kar/Angelman.kt.txt             | Loss of Chr15q11-13                                                               | deletion(Chr15a:2050001-33400000)                                                          |
| ./Precompiled_Kar/Klinefelter.kt.txt          | Additional ChrX                                                                   | chromosomal_duplication(ChrXa)                                                             |
| ./Precompiled_Kar/Marfan_FBN1.kt.txt          | Loss of FBN1 gene region                                                          | deletion(Chr15a:48,408,313-48,645,709)                                                     |
| ./Precompiled_Kar/Marfan_FGFBR2.kt.txt        | Loss of FGFBR2 gene region                                                        | deletion(Chr3:30,606,601-30,694,142)                                                       |
| ./Precompiled_Kar/Hunter_female.kt.txt        | Loss of IDS gene region                                                           | deletion(ChrXa:149,476,988-149,505,306); deletion(ChrXb:149,476,988-149,505,306)           |
| ./Precompiled_Kar/Hunter_male.kt.txt          | Loss of IDS gene region                                                           | deletion(ChrXa:149,476,988-149,505,306)                                                    |
| ./Precompiled_Kar/Achondroplasia_FGFR3.kt.txt | Loss of FGFR3                                                                     | deletion(Chr4a:1,793,293-1,808,867)                                                        |
| ./Precompiled_Kar/Protanopia_female.kt.txt    | Loss of OPN1LW gene region                                                        | deletion(ChrXa:154,144,243-154,159,032); deletion(ChrXb:154,144,243-154,159,032)           |
| ./Precompiled_Kar/Protanopia_male.kt.txt      | Loss of OPN1LW gene region                                                        | deletion(ChrXa:154,144,243-154,159,032)                                                    |
| ./Precompiled_Kar/Deuteranopia_female.kt.txt  | Loss of OPN1MW gene region                                                        | deletion(ChrXa:154,182,596-154,196,861); deletion(ChrXb:154,182,596-154,196,861)           |
| ./Precompiled_Kar/Deuteranopia_male.kt.txt    | Loss of OPN1MW gene region                                                        | deletion(ChrXa:154,182,596-154,196,861)                                                    |


**NOTE1:** template JSON file can be found in `./Sample_JSON/random.json`
- Manual mode: not recommended in later layers as the indices are relative to the current genome, not the 
standard hg38 genome. Manual calculation needs to be done to account for off-shifts.
(e.g. if a deletion is to be made on Chr1:100,000-200,000 on a genome that had a prior deletion
on Chr1:50,000-100,000, the manual parameter needs to address deletion range Chr1:50,000-150,000;
if deletion is made on Chr1:100,000-200,000, then the actual deleted region would be Chr1:50,000-100,000
and Chr1:150,000-250,000): <br />
`python ./Main/KarSimulator.py manual`
#### Arguments:
| Argument | Type | Description                                                                   |
|----------|------|-------------------------------------------------------------------------------|
| `--json` | FILE | JSON file containing a list of SVs to be operated on, to be executed in-order |
**NOTE1:** template JSON file can be found in `./Sample_JSON/manual.json`

**NOTE2:** the required parameters for each SV is the same as addressed in the **supported SV** image

### Output FASTA
`python ./Main/KarSimulator.py fasta`
#### Arguments:
| Argument   | Type | Description                                                                      |
|------------|------|----------------------------------------------------------------------------------|
| `--name`   | STR  | (default: the prefix of the input kar file) <br /> name of the output FASTA file |
| `--genome` | FILE | (default: gh38 in ./Genomes/) Genome FASTA file                                  |
| `--kar`    | FILE | Karyotype file containing the input karyotype                                    |
| `-o`       | FILE | (default: same directory as Karyotype file) output directory                     |
