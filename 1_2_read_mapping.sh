bwa mem -t 4 -R '@RG\tID:sample\tSM:sample\tPL:ILLUMINA' H37Rv.fasta output_forward_paired.fastq output_reverse_paired.fastq > mapped_reads.sam

samtools view -Sb mapped_reads.sam > mapped_reads.bam
samtools sort mapped_reads.bam -o sorted_reads.bam
samtools index sorted_reads.bam

picard MarkDuplicates I=sorted_reads.bam O=dedup_reads.bam M=metrics.txt
