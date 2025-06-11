## Separate a FASTA file into multiple sub-FASTA, each contains the sequence of a single header
## Used for dotplotting in dgenies


input="/home/zhaoyang-new/Genomes/GRCh19/hg19.fa"
output_dir="/home/zhaoyang-new/Genomes/GRCh19/header_splitted/"

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