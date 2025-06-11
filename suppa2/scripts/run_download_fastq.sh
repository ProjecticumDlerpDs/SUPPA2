#!/bin/bash

# Vereist: raw_data/e85_srr_list.txt moet bestaan (zie README voor uitleg)

# Dit script leest de lijst met SRR-nummers en downloadt alle FASTQ bestanden met prefetch + fasterq-dump
# Bestaande bestanden worden overgeslagen om dubbele downloads te voorkomen.

# Bestand met SRR nummers
SRR_LIST="raw_data/e85_srr_list.txt"

# Map waar je FASTQ bestanden in komen
OUTDIR="raw_data/fastq_files"

# Check of het SRR lijstbestand bestaat
if [ ! -f "$SRR_LIST" ]; then
  echo "Fout: bestand met SRR nummers niet gevonden: $SRR_LIST"
  echo "Zorg dat je dit bestand hebt gedownload volgens de instructies in de README."
  exit 1
fi

mkdir -p $OUTDIR

# Loop over elk SRR nummer in je lijst
while read -r SRR; do
    echo "Processing $SRR..."
    
    # Check of de FASTQ bestanden al bestaan (controleer op beide reads als paired-end)
    if [[ -f "$OUTDIR/${SRR}_1.fastq.gz" && -f "$OUTDIR/${SRR}_2.fastq.gz" ]]; then
    echo "$SRR is al gedownload, wordt overgeslagen."
    continue
    fi
    
    # Download het .sra bestand met prefetch
    echo "prefetching $SRR..."
    prefetch "$SRR"
    
    # Zet .sra om aar FASTQ met fasterq-dump
    echo "Converting $SRR to FASTQ..."
    fastq-dump "$SRR" --split-files --gzip --gzip -O "$OUTDIR"
    
    # Verwijder het gedownloade .sra bestand om schijfruimte te besparen
    echo "Verwijderen van .sra bestand voor $SRR..."
    rm -f "$SRR".sra
    
done < "$SRR_LIST"

echo "Download en conversie klaar!"

