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
module load blast/2.10.1

# === USER PATHS ===
BLAST_DIR="/uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/SPAdes_scripts/blast_result/filtered_coi_out"
FASTA_BASE="/uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/SPAdes_scripts/denovo_assembly"
CSV_OUT="consensus_multi_hits.csv"
FASTA_OUT="consensus_multi_hits.fasta"

# === Init output files ===
echo "SampleID,ContigID,HitAccession,PercentIdentity,Evalue,QueryCoverage,ScientificName,ContigSequence" > "$CSV_OUT"
> "$FASTA_OUT"

# === Main loop ===
for FILE in "$BLAST_DIR"/*.filteredcoioutput; do
    SAMPLE=$(basename "$FILE" .filteredcoioutput)
    CONTIG_PATH="$FASTA_BASE/$SAMPLE/contigs.fasta"

    if [[ ! -f "$FILE" ]]; then
        echo "⚠️ No BLAST file for $SAMPLE, skipping."
        continue
    fi

    if [[ ! -f "$CONTIG_PATH" ]]; then
        echo "⚠️ Missing contigs.fasta for $SAMPLE, skipping."
        continue
    fi

    while IFS=$'\t' read -r QSEQID SSEQID PIDENT LENGTH MISMATCH GAPOPEN QSTART QEND SSTART SEND EVALUE BITSCORE QLEN SLEN; do
        QCOV=$(awk "BEGIN { printf \"%.2f\", (${LENGTH} / ${QLEN}) * 100 }")

        # Try to find scientific name from the optional alignment text file (if it exists)
        ALIGN_FILE="${FILE%.filteredcoioutput}-Alignment.txt"
        SCI_NAME="Unknown"
        if [[ -f "$ALIGN_FILE" ]]; then
            SCI_NAME_LINE=$(grep -A 5 "$QSEQID" "$ALIGN_FILE" | tail -1)
            SCI_NAME=$(echo "$SCI_NAME_LINE" | awk '{print $1 " " $2}')
        fi

        # Extract contig sequence from contigs.fasta
        SEQ=$(awk -v id="$QSEQID" '
            BEGIN { found=0; seq="" }
            /^>/ { if (found) exit; found=($0 ~ ">"id) }
            found && !/^>/ { seq = seq $0 }
            END { print seq }
        ' "$CONTIG_PATH")

        if [[ -z "$SEQ" ]]; then
            echo "⚠️ Could not find contig $QSEQID in $SAMPLE"
            continue
        fi

        # Write CSV row
        echo "$SAMPLE,$QSEQID,$SSEQID,$PIDENT,$EVALUE,$QCOV,$SCI_NAME,\"$SEQ\"" >> "$CSV_OUT"

        # Write FASTA entry
        echo ">$SAMPLE|$QSEQID|$SCI_NAME|$SSEQID" >> "$FASTA_OUT"
        echo "$SEQ" >> "$FASTA_OUT"

    done < "$FILE"

done

echo "✅ Done! Output written to:"
echo "   - $CSV_OUT"
echo "   - $FASTA_OUT"

chmod -R g+w "$BLAST_DIR" "$FASTA_BASE" "$CSV_OUT" "$FASTA_OUT"
