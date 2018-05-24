rule process_scores:
    input:
        alignments=expand('alignments/{{run}}-{sample}.bam',
                          sample=REFERNCE_TRANSCRIPTOME_SEQS),
        blacklist='refs/badtranscripts.txt'
    output: 'tables/alignment-scores-{run}.txt'
    run:
        import pandas as pd
        import subprocess as sp
        import numpy as np

        blacklist = set(open(input.blacklist).read().split())
        readids_ambig = set()

        def load_scores(filename):
            nonlocal readids_ambig

            samplename = filename.rsplit('-', 1)[-1].rsplit('.', 1)[0]
            with sp.Popen('samtools view -F4 {} | cut -f1,3,5'.format(filename), shell=True,
                          stdout=sp.PIPE) as stproc:
                tbl = pd.read_table(stproc.stdout, names=['read_id', 'ref', samplename],
                                    dtype={'read_id': str, 'ref': str,
                                           samplename: np.int32})
                readids_ambig |= set(tbl[tbl['ref'].isin(blacklist)]['read_id'].tolist())
                return tbl[['read_id', samplename]].sort_values(by=samplename,
                                    ascending=False).groupby('read_id').first()

        scores = None
        for inpfile in input.alignments:
            print('Loading', inpfile)
            partialtbl = load_scores(inpfile)
            if scores is None:
                scores = partialtbl
            else:
                scores = pd.merge(scores, partialtbl, how='outer',
                                  left_index=True, right_index=True)

        scores = scores.fillna(0).clip(0).astype(np.int32)
        scores = scores[~scores.index.to_series().isin(readids_ambig)].copy()
        best_score = pd.Series(np.sort(scores.values)[:, -1], index=scores.index)
        second_best_score = pd.Series(np.sort(scores.values)[:, -2], index=scores.index)
        best_score_idx = scores.idxmax(axis=1)

        sample2index = {
            refinfo[0]: idxi
            for idxi, refinfo in enumerate(config['reference_transcriptomes'])}
        scores['_best_score'] = best_score
        scores['_best_score_diff'] = best_score - second_best_score
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
            (seqs['_best_score_diff'] >= selection_criteria['min_alignment_score_diff']) &
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

rule subsample_training_inputs:
    input:
        catalog='tables/selected-signal-matches-{run}.txt',
        fulldata='arrays/full-training-{run}.npy'
    output: 'traindata/signals-{run}-s{size}-t{trim}.hdf5'
    run:
        from sklearn.utils import class_weight
        from sklearn.preprocessing import OneHotEncoder
        import h5py
        import numpy as np
        import pandas as pd

        subsamplesize = int(wildcards.size)
        trimsize = int(wildcards.trim)

        fullsig = np.load(input.fulldata)[:, -trimsize:, :]
        allreads = pd.read_table(input.catalog)
        indexcounts = allreads['_best_score_idxi'].value_counts().clip(None, subsamplesize)
        testcounts = (indexcounts *
                      config['train_data_transform']['test_data_split']).astype(np.int64)

        input_signals, input_labels, input_indices = [], [], []
        for idx, cnt in indexcounts.sort_index().items():
            readidx = allreads[allreads['_best_score_idxi'] == idx].index.tolist()
            readidx = np.random.permutation(readidx)[:cnt]
            input_signals.append(fullsig[readidx])
            input_indices.append(readidx)

            sublabels = np.array([idx] * cnt, dtype=np.int16)
            testidx = np.random.permutation(list(range(cnt)))[:testcounts.loc[idx]]
            sublabels[testidx] = -(sublabels[testidx] + 1)
            input_labels.append(sublabels)

        input_signals = np.vstack(input_signals)
        input_labels = np.hstack(input_labels)
        input_indices = np.hstack(input_indices)

        shuforder = np.random.permutation(list(range(len(input_signals))))
        input_signals = input_signals[shuforder]
        input_labels = input_labels[shuforder]
        input_indices = input_indices[shuforder]

        is_test = input_labels < 0
        train_labels = input_labels[~is_test]
        test_labels = -(input_labels[is_test] + 1)

        ohencoder = OneHotEncoder(sparse=False, n_values=len(indexcounts))
        train_onehot = ohencoder.fit_transform(train_labels[:, np.newaxis]).astype(np.float32)
        test_onehot = ohencoder.fit_transform(test_labels[:, np.newaxis]).astype(np.float32)

        clsweight_train = class_weight.compute_class_weight('balanced',
                            list(range(len(indexcounts))), train_labels)
        clsweight_test = class_weight.compute_class_weight('balanced',
                            list(range(len(indexcounts))), test_labels)

        with h5py.File(output[0], 'w') as h5:
            grp = h5.create_group('training')
            grp['signals'] = input_signals[~is_test]
            grp['labels'] = train_labels
            grp['onehot'] = train_onehot
            grp['readid'] = allreads.iloc[input_indices[~is_test]]['read_id'].astype('S36')
            grp['weights'] = clsweight_train

            grp = h5.create_group('testing')
            grp['signals'] = input_signals[is_test]
            grp['labels'] = test_labels
            grp['onehot'] = test_onehot
            grp['readid'] = allreads.iloc[input_indices[is_test]]['read_id'].astype('S36')
            grp['weights'] = clsweight_test


rule train_demux_nn:
    input: 'traindata/signals-{run}-s{size}-t{trim}.hdf5'
    output: 'models/{run}-s{size}-t{trim}-l{preset}/test-predictions.txt'
    threads: 99
    run:
        import json
        import subprocess as sp
        global_params = json.dumps(config['training_parameters']['global'])
        layer_def = json.dumps(config['training_parameters'][wildcards.preset])
        sp.check_call(['scripts/train_demux_nn.py', global_params, layer_def,
                      input[0], os.path.dirname(output[0])])


# ex: syntax=snakemake sw=4 sts=4 et
