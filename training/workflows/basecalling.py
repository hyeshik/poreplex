rule basecall:
    output: '{name}.albacore/.hyeseq_basecall_done'
    threads: 40
    run:
        input_dir = config['data'][wildcards.name]
        output_dir = os.path.dirname(output[0])
        shell('read_fast5_basecaller.py \
                --input "{input_dir}" \
                --resume --worker_threads {threads} \
                --save_path "{output_dir}" \
                --flowcell {config[flowcell]} --kit {config[kit]} \
                --disable_filtering --disable_pings --recursive \
                --files_per_batch_folder 1000 \
                --output_format fastq,fast5 \
                --reads_per_fastq_batch 1000')
        shell('touch {output}')

rule datadir_reorganize:
    input: '{name}.albacore/.hyeseq_basecall_done'
    output: '{name}/.hyeseq_reorganization_done', '{name}/sequences.fasta.gz'
    threads: 16
    shell: 'scripts/albacore_reorganize.py -o {wildcards.name} \
                -i {wildcards.name}.albacore -p {threads} \
            && touch {output}'

rule generate_catalog:
    input: expand('{name}/.hyeseq_reorganization_done', name=RUN_NAMES)
    output: 'sequencing_summary.feather'
    shell: 'scripts/generate_fast5_catalog.py {RUN_NAMES}'
