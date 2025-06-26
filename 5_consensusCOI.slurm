#!/bin/bash
#SBATCH --time=336:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=20          # same as $max set in ForkManager
#SBATCH --account=saarman-np
#SBATCH --partition=saarman-shared-np   
#SBATCH --job-name=consensus-COI
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=emily.calhoun@usu.edu

# Load any needed modules 
module load blast/2.11.0

# Define paths
BLAST_DIR="/uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/SPAdes_scripts/blast_result/filtered_coi_out"                # Directory containing filtered coi 
FASTA_BASE="/uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/SPAdes_scripts/denovo_assembly"                             # Directory containing sample subfolders with contigs.fasta
CSV_OUT="/uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/SPAdes_scripts/blast_result/consensus_multi_hits.csv"
FASTA_OUT="/uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/SPAdes_scripts/blast_result/consensus_multi_hits.fasta"

# === Initialize the output files ===
echo "SampleID,ContigID,HitAccession,PercentIdentity,Evalue,QueryCoverage,ScientificName,ContigSequence" > "$CSV_OUT"
> "$FASTA_OUT"  # Truncate/create the FASTA file

# === Process each BLAST output file ===
for FILE in "$BLAST_DIR"/*.m8; do
    SAMPLE=$(basename "$FILE" .m8)
    CONTIG_PATH="$FASTA_BASE/$SAMPLE/contigs.fasta"

    # Skip if required files don't exist
    if [[ ! -f "$FILE" ]]; then
        echo "⚠️ BLAST file missing for $SAMPLE, skipping."
        continue
    fi
    if [[ ! -f "$CONTIG_PATH" ]]; then
        echo "⚠️ contigs.fasta missing for $SAMPLE, skipping."
        continue
    fi

    ALIGN_FILE="${FILE%.m8}-Alignment.txt"

    # === Process every line (hit) in the BLAST file ===
    while IFS=$'\t' read -r QSEQID SSEQID PIDENT LENGTH MISMATCH GAPOPEN QSTART QEND SSTART SEND EVALUE BITSCORE QLEN SLEN; do
        
        # Calculate query coverage as a percentage
        QCOV=$(awk "BEGIN { printf \"%.2f\", (${LENGTH} / ${QLEN}) * 100 }")

        # Try to get the scientific name from optional alignment summary file
        SCI_NAME="Unknown"
        if [[ -f "$ALIGN_FILE" ]]; then
            SCI_NAME_LINE=$(grep -A 5 "$SSEQID" "$ALIGN_FILE" | tail -1)
            SCI_NAME=$(echo "$SCI_NAME_LINE" | awk '{print $1 " " $2}')
        fi

        # Extract the corresponding contig sequence from contigs.fasta using SSEQID
        SEQ=$(awk -v id="$SSEQID" '
            BEGIN { found=0; seq="" }
            /^>/ { if (found) exit; found=($0 ~ ">"id) }
            found && !/^>/ { seq = seq $0 }
            END { print seq }
        ' "$CONTIG_PATH")

        # Skip if contig not found in FASTA
        if [[ -z "$SEQ" ]]; then
            echo "⚠️ Could not find contig $SSEQID in $SAMPLE"
            continue
        fi

        # === Write to CSV ===
        echo "$SAMPLE,$SSEQID,$QSEQID,$PIDENT,$EVALUE,$QCOV,$SCI_NAME,\"$SEQ\"" >> "$CSV_OUT"

        # === Write to FASTA ===
        echo ">$SAMPLE|$SSEQID|$SCI_NAME|$QSEQID" >> "$FASTA_OUT"
        echo "$SEQ" >> "$FASTA_OUT"

    done < "$FILE"

done

# === Done ===
echo "✅ Complete!"
echo "Output written to:"
echo "  CSV:   $CSV_OUT"
echo "  FASTA: $FASTA_OUT"

chmod -R g+w "$BLAST_DIR" "$FASTA_BASE" "$CSV_OUT" "$FASTA_OUT"
