gatk HaplotypeCaller -R H37Rv.fasta -I dedup_reads.bam -O raw_variants.vcf

gatk VariantFiltration -R H37Rv.fasta -V raw_variants.vcf --filter-expression "QD < 2.0 || MQ < 40.0" --filter-name "LowQual" -O filtered_variants.vcf