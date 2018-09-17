rule preprocess:
    output:
        bookmark_file='{name}/.poreplex_preprocessing_done',
        adapter_inventory='{name}/adapter-dumps/inventory.h5'
    threads: 40
    run:
        input_dir = config['data'][wildcards.name]
        output_dir = os.path.dirname(output.bookmark_file)
        shell('poreplex --input "{input_dir}" --output "{output_dir}" --parallel {threads} \
                --dump-adapter-signals --dump-basecalled-events --trim-adapter \
                --fast5 --symlink-fast5 --basecall')
        shell('touch {output.bookmark_file}')


rule generate_fasta:
    input: '{name}/.poreplex_preprocessing_done'
    output: 'sequences/{name}.fa.gz'
    threads: 2
    shell: """zcat {wildcards.name}/fastq/pass.fastq.gz | \
            awk '{{ if (NR % 4 == 1) print ">" substr($0, 2); \
                 else if (NR % 4 == 2) {{ gsub(/U/, "T"); print $0; }} }}' | \
            gzip -c - > {output}"""

rule generate_catalog:
    input: expand('{name}/.poreplex_preprocessing_done', name=RUN_NAMES)
    output: 'sequencing_summary.feather'
    shell: 'scripts/generate_fast5_catalog.py {RUN_NAMES}'
