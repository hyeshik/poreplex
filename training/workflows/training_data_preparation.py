rule process_scores:
    input: expand('alignments/{{run}}-{sample}.bam', sample=REFERNCE_TRANSCRIPTOME_SEQS)
    output: 'tables/alignment-scores-{run}.txt'
    run:
        import pandas as pd
        import subprocess as sp
        import numpy as np

        def load_scores(filename):
            samplename = filename.rsplit('-', 1)[-1].rsplit('.', 1)[0]
            with sp.Popen('samtools view -F4 {} | cut -f1,4'.format(filename), shell=True,
                          stdout=sp.PIPE) as stproc:
                tbl = pd.read_table(stproc.stdout, names=['read_id', samplename],
                                    dtype={'read_id': str, samplename: np.int32})
                return tbl.sort_values(by=samplename,
                                       ascending=False).groupby('read_id').first()

        scores = None
        for inpfile in input:
            print('Loading', inpfile)
            partialtbl = load_scores(inpfile)
            if scores is None:
                scores = partialtbl
            else:
                scores = pd.merge(scores, partialtbl, how='outer',
                                  left_index=True, right_index=True)

        scores = scores.fillna(0).clip(0).astype(np.int32)
        best_score = pd.Series(np.sort(scores.values)[:, -1], index=scores.index)
        second_best_score = pd.Series(np.sort(scores.values)[:, -2], index=scores.index)
        best_score_idx = scores.idxmax(axis=1)

        sample2index = {
            refinfo[0]: idxi
            for idxi, refinfo in enumerate(config['reference_transcriptomes'])}
        scores['_best_score'] = best_score
        scores['_best_score_ratio'] = best_score / (best_score + second_best_score).clip(1)
        scores['_best_score_idx'] = best_score_idx
        scores['_best_score_idxi'] = best_score_idx.apply(sample2index.__getitem__)

        scores.to_csv(output[0], sep='\t')

rule select_reads:
    input:
        alignment_score='tables/alignment-scores-{run}.txt',
        signal_extract=lambda wc: config['adapter_signal_data'][wc.run]+'/adapter-widths.txt',
        sequencing_summary='sequencing_summary.feather',
    output:
        adapter_signal_matches='tables/adapter-signal-matches-{run}.txt',
        selected_signal_matches='tables/selected-signal-matches-{run}.txt'
    run:
        import pandas as pd
        import feather

        alnscores = pd.read_table(input.alignment_score, index_col=0)
        seqsummary = feather.read_dataframe(input.sequencing_summary)

        fast5prefix = wildcards.run + '/fast5/'

        adapter_signal_list = pd.read_table(input.signal_extract,
                names=['signal_file', 'adapter_start', 'adapter_end'])
        adapter_signal_list['adapter_length'] = (
                adapter_signal_list['adapter_end'] - adapter_signal_list['adapter_start'])
        adapter_signal_list['filename'] = (
                adapter_signal_list['signal_file'].apply(
                    lambda x: fast5prefix + x[:-4] + '.fast5'))

        seqs = (pd.merge(
            alnscores,
            pd.merge(adapter_signal_list, seqsummary, how='inner',
                     left_on='filename', right_on='filename'),
            how='inner', left_index=True, right_on='read_id')
                .sort_values(by=['_best_score_idxi', 'read_id']))
        seqs.to_csv(output.adapter_signal_matches, sep='\t', index=False)

        selection_criteria = config['train_data_selection']
        selected = seqs[
            (seqs['_best_score'] >= selection_criteria['min_alignment_score']) &
            (seqs['_best_score_ratio'] >= selection_criteria['min_alignment_score_ratio']) &
            (seqs['adapter_length'] >= selection_criteria['min_adapter_length']) &
            (seqs['adapter_length'] <= selection_criteria['max_adapter_length']) &
            (seqs['sequence_length_template'] >= selection_criteria['min_sequence_length'])]

        selected.to_csv(output.selected_signal_matches, sep='\t', index=False)

rule prepare_training_data:
    input:
        catalog='tables/selected-signal-matches-{run}.txt',
        signal_extract=lambda wc: config['adapter_signal_data'][wc.run]+'/adapter-widths.txt'
    output: 'arrays/full-training-{run}.npy'
    threads: 8
    run:
        import os
        signal_dir = os.path.dirname(input.signal_extract)
        cfg = config['train_data_transform']
        shell('scripts/prepare_training_data.py {signal_dir} \
                {cfg[trim_adapter_length]} {cfg[summarization_window_size]} \
                {input.catalog} {output} {threads}')

