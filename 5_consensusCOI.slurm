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

# Paths
FILTERED_DIR="/uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/SPAdes_scripts/blast_result/filtered_coi_out"
ASSEMBLY_DIR="/uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/SPAdes_scripts/denovo_assembly"
CSV_OUT="/uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/SPAdes_scripts/blast_result/consensus_top_hits.csv"
FASTA_OUT="/uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/SPAdes_scripts/blast_result/consensus_top_hits.fasta"

# Init output files
echo "SampleID,ContigID,TopHitAccession,PercentIdentity,Evalue,QueryCoverage,ScientificName,ContigSequence" > "$CSV_OUT"
> "$FASTA_OUT"

# Loop through filtered BLAST results
for FILE in "$FILTERED_DIR"/*.filteredcoioutput; do
    SAMPLEID=$(basename "$FILE" .filteredcoioutput)
    TOPHIT=$(head -n 1 "$FILE")

    if [[ -z "$TOPHIT" ]]; then
        echo "Skipping $SAMPLEID: no top hit"
        continue
    fi

    # Parse top hit
    CONTIGID=$(echo "$TOPHIT" | cut -f1)
    ACCESSION=$(echo "$TOPHIT" | cut -f2)
    PID=$(echo "$TOPHIT" | cut -f3)
    ALIGNLEN=$(echo "$TOPHIT" | cut -f4)
    EVAL=$(echo "$TOPHIT" | cut -f11)
    QLEN=$(echo "$TOPHIT" | cut -f13)
    QCOV=$(awk -v a="$ALIGNLEN" -v q="$QLEN" 'BEGIN {printf "%.2f", (a/q)*100}')

    # Get scientific name from *_Alignment.txt
    SCI_NAME=$(grep "$ACCESSION" "$FILTERED_DIR"/../*_Alignment.txt | tail -1 | awk '{print $1 " " $2}')
    [[ -z "$SCI_NAME" ]] && SCI_NAME="Unknown"

    # Get contig path
    CONTIG_PATH="$ASSEMBLY_DIR/$SAMPLEID/contigs.fasta"
    if [[ ! -f "$CONTIG_PATH" ]]; then
        echo "Missing contigs.fasta for $SAMPLEID"
        continue
    fi

    # Extract contig sequence
    SEQ=$(awk -v id=">${CONTIGID}" '
        $0 ~ id {flag=1; next}
        /^>/ {flag=0}
        flag {seq=seq $0}
        END {print seq}' "$CONTIG_PATH")

    if [[ -z "$SEQ" ]]; then
        echo "No sequence found for $CONTIGID in $SAMPLEID"
        continue
    fi

    # Write to CSV (quoted sequence to preserve line integrity)
    echo "$SAMPLEID,$CONTIGID,$ACCESSION,$PID,$EVAL,$QCOV,$SCI_NAME,\"$SEQ\"" >> "$CSV_OUT"

    # Write to multi-FASTA
    echo ">$SAMPLEID|$CONTIGID|$SCI_NAME" >> "$FASTA_OUT"
    echo "$SEQ" >> "$FASTA_OUT"
done

echo "✅ All done. Output written to:"
echo "- $CSV_OUT"
echo "- $FASTA_OUT"
