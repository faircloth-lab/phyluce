# Check assemblies

## Spades

```bash
python ../bin/assembly/phyluce_assembly_assemblo_spades \
	--config ../phyluce/tests/test-conf/assembly.conf \
	--cores 3 --output spades
```


## ABYSS

* we need to force abyss-se because regular abyss-pe will fail:

```bash
python ../bin/assembly/phyluce_assembly_assemblo_abyss \
	--config ../phyluce/tests/test-conf/assembly.conf \
	--cores 3 \
	--output abyss \
	--abyss-se
```

## Velvet

```bash
python ../bin/assembly/phyluce_assembly_assemblo_velvet \
	--config ../phyluce/tests/test-conf/assembly.conf \
	--cores 3 \
	--output velvet
```

# Find UCE contigs

```bash
python ../bin/assembly/phyluce_assembly_match_contigs_to_probes 
	--contigs spades/contigs \
	--probes ../phyluce/tests/probes/uce-5k-probes.fasta \
	--output probe-match
```

* expecting:

```
alligator_mississippiensis: 24
gallus_gallus: 20
peromyscus_maniculatus: 16
rana_sphenocephafa: 6
```

# Get match counts

```bash
python ../bin/assembly/phyluce_assembly_get_match_counts \
	--locus-db probe-match/probe.matches.sqlite \
	--taxon-list-config ../phyluce/tests/test-conf/taxon-set.conf \
	--taxon-group all \
	--incomplete-matrix \
	--output taxon-set.conf
```

* expecting:

```
There are 25 UCE loci in an INCOMPLETE matrix
```

# Get FASTAs

```bash
python ../bin/assembly/phyluce_assembly_get_fastas_from_match_counts \
	--contigs spades/contigs \
	--locus-db probe-match/probe.matches.sqlite \
	--match-count-output taxon-set.conf \
	--incomplete-matrix taxon-set.incomplete \
	--output taxon-set.fasta
```

* expecting:

```
	There are 24 UCE loci for alligator_mississippiensis
	There are 20 UCE loci for gallus_gallus
	There are 16 UCE loci for peromyscus_maniculatus
	There are 6 UCE loci for rana_sphenocephafa
```

# Align checks

* both produce 15 alignments:

```bash
	python ../bin/align/phyluce_align_seqcap_align \
	    --fasta taxon-set.fasta \
	    --output mafft \
	    --taxa 4 \
	    --aligner mafft \
	    --cores 4 \
	    --incomplete-matrix \
	    --output-format fasta \
	    --no-trim

	python ../bin/align/phyluce_align_seqcap_align \
	    --fasta taxon-set.fasta \
	    --output muscle \
	    --taxa 4 \
	    --aligner muscle \
	    --cores 4 \
	    --incomplete-matrix \
	    --output-format fasta \
	    --no-trim
```


# Check trimming

```bash
python ../bin/align/phyluce_align_get_gblocks_trimmed_alignments_from_untrimmed \
    --alignments mafft \
    --output mafft-gblocks \
    --cores 4

python ../bin/align/phyluce_align_get_trimal_trimmed_alignments_from_untrimmed \
    --alignments mafft \
    --output mafft-trimal \
    --cores 4

python ../bin/align/phyluce_align_get_gblocks_trimmed_alignments_from_untrimmed \
    --alignments muscle \
    --output muscle-gblocks \
    --cores 4

python ../bin/align/phyluce_align_get_trimal_trimmed_alignments_from_untrimmed \
    --alignments muscle \
    --output muscle-trimal \
    --cores 4
```

# Remove locus Names

```bash
python ../bin/align/phyluce_align_remove_locus_name_from_nexus_lines \
    --alignments mafft-trimal \
    --output mafft-trimal-clean \
    --cores 4

python ../bin/align/phyluce_align_remove_locus_name_from_nexus_lines \
    --alignments mafft-gblocks \
    --output mafft-gblocks-clean \
    --cores 4
```

# Get loci with min taxa

```bash
python ../bin/align/phyluce_align_get_only_loci_with_min_taxa \
    --alignments mafft-gblocks-clean \
    --taxa 4 \
    --percent 0.75 \
    --output mafft-gblocks-clean-75p \
    --cores 4
```

expected:

```
Copied 15 alignments of 15 total containing ≥ 0.75 proportion of taxa (n = 3)
```

# Prep for RAXML

```bash
python ../bin/align/phyluce_align_format_nexus_files_for_raxml \
    --alignments mafft-gblocks-clean-75p \
    --output mafft-gblocks-clean-75p-raxml \
    --charsets
```

* expected phylip header:

```
4 9257
```