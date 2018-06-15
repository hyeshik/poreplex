rule make_minimap2_index:
    input: 'refs/{name}.fa.gz'
    output: 'refs/{name}.mmidx'
    shell: 'zcat {input} | minimap2 -H -k 13 -d {output} /dev/stdin'

rule map_to_controls:
    input: seq='sequences/{name}.fa.gz', idx='refs/controls.mmidx'
    output: 'alignments/{name}-controls.bam'
    threads: 40
    shell: 'minimap2 -a --cs=long -t {threads} -k 13 -H {input.idx} {input.seq} | \
            samtools view -@ 8 -b -o {output} /dev/stdin'

rule map_to_transcriptome:
    input:
        seq='sequences/{name}.fa.gz',
        idx=lambda wc: REFERNCE_TRANSCRIPTOME_SEQS[wc.refsample].replace('.fa.gz', '.mmidx')
    output:
        'alignments/{{name}}-{{refsample,{0}}}.bam'.format('|'.join(REFERNCE_TRANSCRIPTOME_SEQS))
    threads: 40
    shell: 'minimap2 -a --cs=long -t {threads} -k 13 -H {input.idx} {input.seq} | \
            samtools view -@ 8 -b /dev/stdin | \
            samtools sort -@ {threads} -o {output} /dev/stdin'

rule make_reference_mapq_table:
    input: 'alignments/{name}-controls.bam'
    output: 'tables/controls-mapq-{name}.txt.gz'
    threads: 5
    shell: '(echo "mapq\treadid"; \
             samtools view -F4 {input} | cut -f1,5 | sort -k1,1 -k2,2rn | \
             awk -F\'\t\' \'BEGIN {{ OFS="\\t"; }} {{ print $2, $1; }}\' | \
             uniq -f1) | gzip -c - > {output}'

rule index_bam:
    input: 'alignments/{name}.bam'
    output: 'alignments/{name}.bam.bai'
    shell: 'samtools index {input}'
