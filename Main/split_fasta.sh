## Separate a FASTA file into multiple sub-FASTA, each contains the sequence of a single header
## Used for dotplotting in dgenies


input="/media/zhaoyang-new/workspace/Molecular_Karyotype/KarSimulator/DEMO-Terminal_SV_Simulation/hg19_demo_files/hg19_46XX_demo1.fasta"
output_dir="/media/zhaoyang-new/workspace/Molecular_Karyotype/KarSimulator/DEMO-Terminal_SV_Simulation/hg19_demo_files/demo1_header_split/"

mkdir -p "$output_dir"  # create output directory if it doesn't exist

awk -v outdir="$output_dir" '
    /^>/ {
        if (out) close(out)
        out = substr($0, 2)
        gsub(/[^A-Za-z0-9_.-]/, "_", out)  # sanitize filename
        out = outdir "/" out ".fasta"
    }
    { print >> out }
' "$input"