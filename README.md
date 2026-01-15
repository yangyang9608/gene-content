# gene-content


``
conda create -n agat_env -c conda-forge -c bioconda agat
conda activate agat_env
```

Input

A directory containing GFF files, one sample per file, e.g.

``
Aara.gff
Achr.gff
Acol.gff
...

``

Optional genome FASTA (only needed if you enable phase fixing):

Same directory as GFFs or provided via --genome_dir

Must share the same prefix as the GFF:

Aara.gff → Aara.fa (or .fasta, .fna, .fa.gz, etc.)


#Output structure#

For each sample <sample>, outputs are placed in:

```
out/<sample>/
  run.log
  00.norm.gff3
  01.fix_cds_overlap.gff3
  02.longest_isoform.gff3
  03.phase_fixed.gff3   (only if phase fixing enabled and genome found)
  final.gff3

```
Global outputs:

```
out/summary_agat.tsv
out/failed.list

```
Usage
Basic run (AGAT only)
```
./00.preQC_AGAT.sh -i . -p "A*.gff" -o out

```
Resume behavior (skip existing outputs)

This pipeline is designed to be re-run safely:

If out/<sample>/00.norm.gff3 already exists and is non-empty, the “convert” step is skipped

Same for overlap-fix, keep-longest, and optional phase-fix outputs


