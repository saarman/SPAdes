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

# Define paths
FILTERED_DIR="/uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/SPAdes_scripts/blast_result/filtered_coi_out"
ASSEMBLY_DIR="/uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/SPAdes_scripts/denovo_assembly"
OUTPUT_FILE="consensus_top_hits_with_sequences.csv"

# Write header
echo "SampleID,ContigID,TopHitAccession,PercentIdentity,Evalue,QueryCoverage,ScientificName,ContigSequence" > "$OUTPUT_FILE"

# Loop through filtered BLAST result files
for file in "$FILTERED_DIR"/*.filteredcoioutput; do
    FILENAME=$(basename "$file")
    SAMPLEID="${FILENAME%%.filteredcoioutput}"

    # Read top hit
    TOPHIT=$(head -n 1 "$file")

    if [[ -z "$TOPHIT" ]]; then
        echo "Skipping $SAMPLEID (no top hit)"
        continue
    fi

    # Parse BLAST fields
    CONTIGID=$(echo "$TOPHIT" | cut -f1)
    ACCESSION=$(echo "$TOPHIT" | cut -f2)
    PID=$(echo "$TOPHIT" | cut -f3)
    ALIGNLEN=$(echo "$TOPHIT" | cut -f4)
    EVAL=$(echo "$TOPHIT" | cut -f11)
    QLEN=$(echo "$TOPHIT" | cut -f13)
    QCOV=$(awk -v len="$ALIGNLEN" -v qlen="$QLEN" 'BEGIN {printf "%.2f", (len/qlen)*100}')

    # Get scientific name
    SCI_NAME=$(grep "$ACCESSION" "$FILTERED_DIR"/../*_Alignment.txt | tail -1 | awk '{print $1 " " $2}')
    [[ -z "$SCI_NAME" ]] && SCI_NAME="Unknown"

    # Locate contigs.fasta and extract matching contig sequence
    CONTIG_PATH="$ASSEMBLY_DIR/$SAMPLEID/contigs.fasta"
    if [[ ! -f "$CONTIG_PATH" ]]; then
        echo "Contigs file missing for $SAMPLEID"
        continue
    fi

    # Extract FASTA sequence for matching contig
    SEQUENCE=$(awk -v id=">${CONTIGID}" '
        $0 ~ id {flag=1; next} 
        /^>/ {flag=0} 
        flag {seq=seq $0} 
        END {print seq}' "$CONTIG_PATH")

    [[ -z "$SEQUENCE" ]] && SEQUENCE="NotFound"

    # Output CSV line
    echo "$SAMPLEID,$CONTIGID,$ACCESSION,$PID,$EVAL,$QCOV,$SCI_NAME,\"$SEQUENCE\"" >> "$OUTPUT_FILE"
done

echo "✅ Done! FASTA-inclusive file written to: $OUTPUT_FILE"
