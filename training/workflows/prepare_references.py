rule align_trome2trome:
    input:
        query='refs/{query}-all.fa.gz',
        target='refs/{target}-all.mmidx'
    output: 'refs-filtering/{query}-{target}.bam'
    threads: 40
    run:
        if wildcards.query == wildcards.target:
            shell('touch {output}')
        else:
            shell('minimap2 -a -t {threads} -k 13 -H {input.target} {input.query} | \
                   samtools view -@ 8 -b -o {output} /dev/stdin')

rule make_transcript_blacklist:
    input: 'refs-filtering/{query}-{target}.bam'
    output: 'refs-filtering/badtranscripts-{query}-{target}.txt'
    run:
        import subprocess as sp
        import re

        pat_cigar = re.compile(r'(\d+)M')
        def count_matches(cigar):
            matches = pat_cigar.findall(cigar)
            return sum(map(int, matches)) if matches else 0

        if wildcards.query == wildcards.target:
            shell('touch {output}')
        else:
            with sp.Popen(['samtools', 'view', '-F4', input[0]], stdout=sp.PIPE) as proc, \
                        open(output[0], 'w') as out:
                for line in proc.stdout:
                    fields = line.decode().split('\t')
                    print(fields[0], fields[2], count_matches(fields[5]), sep='\t', file=out)

rule finalize_transcript_blacklist:
    input:
        aligned=expand('refs-filtering/badtranscripts-{query}-{target}.txt',
                       query=REFERNCE_TRANSCRIPTOME_SEQS,
                       target=REFERNCE_TRANSCRIPTOME_SEQS),
        curated='refs/curated-blacklist.txt'
    output: 'refs/badtranscripts.txt'
    run:
        import subprocess as sp
        import pandas as pd
        from io import StringIO

        minlength = (config['reference_transcriptome_filtering']
                           ['blacklist_alignment_match_length'])
        allcontent = StringIO(sp.check_output(['cat'] + input.aligned).decode())
        matches = pd.read_table(allcontent, names=['query', 'target', 'm_length'])
        matches = matches[matches['m_length'] >= minlength]

        foundbadtr = set(matches['query'].tolist() + matches['target'].tolist())
        curatedbadtr = open(input.curated).read().split()
        if curatedbadtr == ['']:
            curatedbadtr.pop()
        print(*sorted(foundbadtr | set(curatedbadtr)), sep='\n', file=open(output[0], 'w'))

