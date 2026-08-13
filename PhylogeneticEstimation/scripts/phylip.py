from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord

alignment = AlignIO.read("temp-align.fa", "fasta")

shortened_records = []
for record in alignment:
    short_name = record.id[:10]
    new_record = SeqRecord(record.seq, id=short_name, description="")
    shortened_records.append(new_record)

shortened_alignment = MultipleSeqAlignment(shortened_records)

AlignIO.write(shortened_alignment, "bam-align.phy", "phylip")
