rule train_demux_nn:
    input: 'traindata/signals-{run}-s{size}-t{trim}.hdf5'
    output:
        testpreds='models/{run}-s{size}-t{trim}-l{preset}/test-predictions.txt',
        model='models/{run}-s{size}-t{trim}-l{preset}/final-model.hdf5'
    threads: 99
    run:
        import json
        import subprocess as sp
        global_params = json.dumps(config['training_parameters']['global'])
        sp.check_call(['scripts/train_demux_nn.py', global_params, wildcards.preset,
                      input[0], os.path.dirname(output[0])])

rule compute_calibration_table:
    input: 'models/{setname}/test-predictions.txt'
    output: 'models/{setname}/score-calibration.txt'
    shell: 'python3 scripts/compute_score_calibration_table.py {input} {output}'

rule put_model_params:
    input:
        model='models/{setname}/final-model.hdf5',
        calibration_table='models/{setname}/score-calibration.txt'
    output: 'models/{setname}/tagged-model.hdf5'
    run:
        shell('cp {input.model} {output}')

        import h5py
        import pandas
        import numpy

        crosscont_penalty = (config['training_parameters']['global']
                                   ['metrics']['cross_contamination_penalty'])
        num_decoy_classes = 1
        num_classes = num_decoy_classes + len(config['reference_transcriptomes'])
        loss_weights = numpy.ones([num_classes, num_classes], dtype=numpy.float32)
        loss_weights[num_decoy_classes:, num_decoy_classes:] *= crosscont_penalty

        caltbl = pandas.read_table(input.calibration_table)
        with h5py.File(output[0], 'r+') as mf:
            grp = mf.create_group('poreplex_params')
            caltblentry = grp.create_dataset('calibration', (len(caltbl),),
                                             dtype=[('phred', '<i4'),
                                                    ('pred_score', '<f8')])
            caltblentry['phred'] = caltbl['phred']
            caltblentry['pred_score'] = caltbl['pred_score']
            grp['loss_weights'] = loss_weights

# ex: syntax=snakemake sw=4 sts=4 et
