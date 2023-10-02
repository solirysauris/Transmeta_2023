import re
from Bio import SeqIO

fasta_file = "/Users/margaritatis/Documents/itmo/meta/virus1.fa"
sequence = SeqIO.parse(fasta_file, "fasta")

cds_list = []
seen_positions = set()

for seq_record in sequence:
    sequence = seq_record.seq
    start_codon = "ATG"
    stop_codons = ["TAA", "TAG", "TGA"]
    cds_pattern = re.compile(rf"{start_codon}.*?(?:{'|'.join(stop_codons)})")

    matches = cds_pattern.findall(str(sequence))

    for match in matches:
        start_pos = sequence.find(match) + 1
        end_pos = start_pos + len(match) - 1

        if (start_pos, end_pos) in seen_positions:
            continue

        chromosome = seq_record.id

        cds_sequence = sequence[start_pos - 1:end_pos]
        rna_sequence = cds_sequence.transcribe()
        amino_acid_sequence = rna_sequence.translate()

        bed_line = f"{chromosome}\t{start_pos}\t{end_pos}\t{cds_sequence}\t{rna_sequence}\t{amino_acid_sequence}\n"

        cds_list.append(bed_line)

        seen_positions.add((start_pos, end_pos))

bed_file = "/Users/margaritatis/Documents/itmo/meta/output.bed"
with open(bed_file, "w") as file:
    file.writelines(cds_list)
