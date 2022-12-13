import os
import tempfile
from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=False, stderr=True)

with tempfile.TemporaryDirectory() as tmpdir:
    shell(
        "STARlong --genomeDir {snakemake.input.genome} "
        "--readFilesIn {snakemake.input.data} "
        "--outFileNamePrefix {tmpdir}/ "
        "--outTmpDir {tmpdir}/STARtmp "
        "--runThreadN {snakemake.threads} "
        "--outFilterMultimapScoreRange 20 "
        "--outFilterScoreMinOverLread 0 "
        "--outFilterMatchNminOverLread 0.66 "
        "--outFilterMismatchNmax 1000 "
        "--winAnchorMultimapNmax 200 "
        "--seedSearchStartLmax 12 "
        "--seedPerReadNmax 100000 "
        "--seedPerWindowNmax 100 "
        "--alignTranscriptsPerReadNmax 100000 "
        "--alignTranscriptsPerWindowNmax 10000 "
        "--outStd BAM_SortedByCoordinate "
        "{snakemake.params.zcat} > {snakemake.output}; "
    )
    
    shell("cat {tmpdir}/Log.out > {snakemake.log}")

