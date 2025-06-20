## Example usage to create karyotypes with terminal tandem duplications and duplication inversions

Intended Goal:
- Create terminal tandem duplications and duplication inversions
- Generate SVs only on one haploid genome
- Left duplication inversion only on p-terminal and right duplication inversion only on q-terminal
- No compounded events

Preparation:
- Download hg19 genome and unzip it
- Recommend downloading from https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/ `hg19.fa.gz` to avoid chromosome name issue

Please run all code exactly at the `Karsimulator/` directory.

Steps
1. (Optional) Initialize WT haploid genome template. 23X and 23Y are already prepared with the following commands. Optionally, you can create haploid genome or genomes with specific set of chromosomes. (e.g. when you want to test functionality on a specific chromosome only). 
   1. `python KarSimulator.py rawGenome --name hg19_WT_23X --sex [ChrX] --copy 1 --index Genomes/hg19_index.txt -o DEMO-Terminal_SV_Simulation/`
   2. `python KarSimulator.py rawGenome --name hg19_WT_23Y --sex [ChrY] --copy 1 --index Genomes/hg19_index.txt -o DEMO-Terminal_SV_Simulation/`
2. (Optional) Modify the `Terminal_SV_configs.JSON` with the desired random parameters. Current setup is the following:
   1. tandem-duplications and duplication-inversions ONLY, with 50-50 probability
   2. all SVs are terminal SVs
   3. a total of 20 SVs per karyotype
   4. duplication-inversions are left foldback when on p-terminal and right foldback when on q-terminal
   5. sizes of event are between 200,000 and 5,000,000 bp, and it is sampled from a uniform distribution
   6. no masking region selected (if given a masking region file, no SV will intersect with regions from this file)
   7. no compounding events
3. Run random SV generator by `python KarSimulator.py random --json DEMO-Terminal_SV_Simulation/Terminal_SV_configs.JSON`
4. Select desired karyotype to output the FASTA by `python KarSimulator.py fasta --name <output file name prefix> --genome <path to hg19 .fa/.fasta genome> --kar <path to karyotype file> -o <output DIR>`
   1. e.g. `python Karsimulator.py fasta --name demo1 --genome Genomes/hg19.fa --kar DEMO-Terminal_SV_Simulation/hg19_demo_files/terminal_duplications_46XX_r1.kt.txt -o DEMO-Terminal_SV_Simulation/hg19_demo_files/`
5. Generate the WT haploid 23X FASTA, run your read simulator separately on the two haploids and combine the coverage

## Optional Verification
You can verify each FASTA by installing and running Dgenies to align and dotplot each chromosome against the reference assembly. This was already done during testing to ensure each SV is indeed present in the FASTA file as intended (https://dgenies.toulouse.inra.fr/). If the fasta files are too big, you can use `Main/split_fasta.sh` to split each fasta into separate chromosomes for alignment.