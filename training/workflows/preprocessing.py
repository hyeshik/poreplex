rule preprocess:
    output: '{name}/.octopus_preprocessing_done'
    threads: 40
    run:
        input_dir = config['data'][wildcards.name]
        output_dir = os.path.dirname(output[0])
        shell('octopus --input "{input_dir}" --output "{output_dir}" --parallel {threads} \
                --dump-adapter-signals --dump-basecalled-events --trim-adapter \
                --fast5 --symlink-fast5 --albacore-onthefly')
        shell('touch {output}')


rule generate_fasta:
    input: '{name}/.octopus_preprocessing_done'
    output: 'sequences/{name}.fa.gz'
    threads: 2
    shell: """zcat {wildcards.name}/fastq/pass.fastq.gz | \
            awk '{{ if (NR % 4 == 1) print ">" substr($0, 2); \
                 else if (NR % 4 == 2) {{ gsub(/U/, "T"); print $0; }} }}' | \
            gzip -c - > {output}"""

rule generate_catalog:
    input: expand('{name}/.octopus_preprocessing_done', name=RUN_NAMES)
    output: 'sequencing_summary.feather'
    shell: 'scripts/generate_fast5_catalog.py {RUN_NAMES}'
