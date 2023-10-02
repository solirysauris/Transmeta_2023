from Bio import SeqIO

fasta_file = "/Users/margaritatis/Documents/itmo/meta/human_T1.fa"
sequence = SeqIO.parse(fasta_file, "fasta")

import re

cds_list = []
seen_positions = set()

for seq_record in sequence:
    sequence = str(seq_record.seq)
    start_codon = "ATG"
    stop_codons = ["TAA", "TAG", "TGA"]
    cds_pattern = re.compile(rf"{start_codon}.*?(?:{'|'.join(stop_codons)})")

    matches = cds_pattern.findall(sequence)

    for match in matches:
        start_pos = sequence.find(match) + 1
        end_pos = start_pos + len(match) - 1

        if (start_pos, end_pos) in seen_positions:
            continue

        chromosome = seq_record.id

        cds_sequence = match

        bed_line = f"{chromosome}\t{start_pos}\t{end_pos}\t{cds_sequence}\n"
        cds_list.append(bed_line)

        seen_positions.add((start_pos, end_pos))

bed_file = "/Users/margaritatis/Documents/itmo/meta/human1.bed"
with open(bed_file, "w") as file:
    file.writelines(cds_list)