# Check assemblies

## Spades

```bash
python ../bin/assembly/phyluce_assembly_assemblo_spades \
	--config ../phyluce/tests/test-conf/assembly.conf \
	--cores 3 --output spades
```

* explode:

```bash
for i in spades/contigs/*.fasta;
do
    ../bin/assembly/phyluce_assembly_get_fasta_lengths --input $i --csv;
done
```

* expected:

```
alligator_mississippiensis.contigs.fasta,24,20187,841.125,10.3267824899,765,1016,839.0,1
gallus_gallus.contigs.fasta,20,17900,895.0,25.6631006864,691,1151,895.5,3
peromyscus_maniculatus.contigs.fasta,21,10280,489.523809524,40.3433253491,210,824,551.0,0
rana_sphenocephafa.contigs.fasta,6,2826,471.0,48.755170666,293,632,473.0,0
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

* explode:

```bash
for i in abyss/contigs/*.fasta;
do
    ../bin/assembly/phyluce_assembly_get_fasta_lengths --input $i --csv;
done
```

* expected:

```
alligator_mississippiensis.contigs.fasta,24,7282,303.416666667,25.4645901251,199,722,259.5,0
gallus_gallus.contigs.fasta,20,8329,416.45,38.1168200888,162,718,410.0,0
peromyscus_maniculatus.contigs.fasta,15,4719,314.6,14.6839076154,230,426,305.0,0
rana_sphenocephafa.contigs.fasta,7,2544,363.428571429,61.7320851757,101,578,351.0,0
```

## Velvet

```bash
python ../bin/assembly/phyluce_assembly_assemblo_velvet \
	--config ../phyluce/tests/test-conf/assembly.conf \
	--cores 3 \
	--output velvet
```

* explode:

```bash
for i in velvet/contigs/*.fasta;
do
    ../bin/assembly/phyluce_assembly_get_fasta_lengths --input $i --csv;
done
```

* expected:

```
alligator_mississippiensis.contigs.fasta,24,18458,769.083333333,19.5311824878,544,997,802.0,0
gallus_gallus.contigs.fasta,22,15102,686.454545455,40.4290157297,237,995,725.0,0
peromyscus_maniculatus.contigs.fasta,16,8392,524.5,29.6086980464,282,683,551.5,0
rana_sphenocephafa.contigs.fasta,7,2645,377.857142857,74.3676288293,103,632,374.0,0
```

# Find UCE contigs

```bash
python ../bin/assembly/phyluce_assembly_match_contigs_to_probes 
	--contigs spades/contigs \
	--probes ../phyluce/tests/probes/uce-5k-probes.fasta \
	--output probe-match
```

* expected

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

* expected

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

## mafft

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
```

* get stats:

```bash
python ../bin/align/phyluce_align_get_align_summary_data \
	--alignments mafft \
	--input-format fasta \
	--output-stats mafft.stats
```

* expected

```
aln,length,sites,differences,characters,gc content,gaps,a count, c count, g count, t count
uce-7014.fasta,1010,17,320,831,40.65,1430,729,548,513,820
uce-7542.fasta,962,0,159,832,33.39,505,761,343,452,825
uce-6095.fasta,932,0,139,840,40.86,527,625,473,454,717
uce-3046.fasta,946,18,228,841,36.05,838,852,532,530,1032
uce-553.fasta,1221,3,217,1011,39.73,1365,1055,727,671,1066
uce-1484.fasta,841,0,126,784,36.82,432,664,419,351,657
uce-865.fasta,903,0,106,810,39.43,386,650,464,452,757
uce-2120.fasta,988,17,272,830,38.56,1027,904,532,596,893
uce-4179.fasta,907,43,444,795,46.03,880,823,618,647,660
uce-6187.fasta,1000,0,115,764,41.62,696,722,498,461,623
uce-6550.fasta,938,0,182,737,36.26,608,738,394,406,668
uce-2539.fasta,850,0,91,742,40.92,407,496,447,430,770
uce-2349.fasta,1026,0,160,799,36.69,677,699,429,452,821
uce-1985.fasta,1013,0,179,819,38.39,606,737,468,466,762
uce-1732.fasta,1049,4,118,852,40.26,1270,950,470,708,798
```

## muscle

```bash
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

* get stats:

```bash
python ../bin/align/phyluce_align_get_align_summary_data \
	--alignments muscle \
	--input-format fasta \
	--output-stats muscle.stats
```

* expected

```
aln,length,sites,differences,characters,gc content,gaps,a count, c count, g count, t count
uce-7014.fasta,1043,15,315,793,40.65,1562,729,548,513,820
uce-7542.fasta,976,0,134,819,33.39,547,761,343,452,825
uce-6095.fasta,932,0,138,840,40.86,527,625,473,454,717
uce-3046.fasta,952,21,220,833,36.05,862,852,532,530,1032
uce-553.fasta,1207,3,215,1013,39.73,1309,1055,727,671,1066
uce-1484.fasta,841,0,119,782,36.82,432,664,419,351,657
uce-865.fasta,903,0,100,810,39.43,386,650,464,452,757
uce-2120.fasta,973,15,268,820,38.56,967,904,532,596,893
uce-4179.fasta,925,36,396,827,46.03,952,823,618,647,660
uce-6187.fasta,1000,0,112,764,41.62,696,722,498,461,623
uce-6550.fasta,926,0,179,740,36.26,572,738,394,406,668
uce-2539.fasta,850,0,87,742,40.92,407,496,447,430,770
uce-2349.fasta,1028,0,154,799,36.69,683,699,429,452,821
uce-1985.fasta,1014,0,150,810,38.39,609,737,468,466,762
uce-1732.fasta,1049,4,116,852,40.26,1270,950,470,708,798
```

# Check trimming

## mafft-gblocks

```bash
python ../bin/align/phyluce_align_get_gblocks_trimmed_alignments_from_untrimmed \
    --alignments mafft \
    --output mafft-gblocks \
    --cores 4
```

* get stats:

```bash
python ../bin/align/phyluce_align_get_align_summary_data \
	--alignments mafft-gblocks \
	--output-stats mafft.gblocks.stats
```

* expected

```
aln,length,sites,differences,characters,gc content,gaps,a count, c count, g count, t count
uce-1732.nexus,609,4,93,609,41.44,194,723,358,571,590
uce-2349.nexus,652,0,131,652,35.64,76,558,328,342,652
uce-1985.nexus,716,0,161,716,37.83,118,601,395,373,661
uce-6550.nexus,679,0,157,679,37.38,167,599,338,361,572
uce-2539.nexus,557,0,79,557,40.36,6,382,331,341,611
uce-6187.nexus,591,0,105,591,41.58,51,518,359,357,488
uce-865.nexus,631,0,78,631,39.79,23,510,380,364,616
uce-2120.nexus,598,17,192,598,38.82,151,675,405,465,696
uce-4179.nexus,569,37,296,569,45.03,213,617,462,467,517
uce-1484.nexus,459,0,59,459,35.84,7,438,279,212,441
uce-7542.nexus,796,0,146,796,33.5,209,690,317,413,759
uce-6095.nexus,510,0,83,510,42.76,17,384,343,304,482
uce-3046.nexus,678,16,181,678,34.93,233,741,435,431,872
uce-553.nexus,772,3,150,772,41.6,273,838,588,583,806
uce-7014.nexus,440,17,128,440,43.27,163,422,370,321,484
```

## mafft-trimal

```bash
python ../bin/align/phyluce_align_get_trimal_trimmed_alignments_from_untrimmed \
    --alignments mafft \
    --output mafft-trimal \
    --cores 4
```

* get stats:

```bash
python ../bin/align/phyluce_align_get_align_summary_data \
	--alignments mafft-trimal \
	--output-stats mafft.trimal.stats
```

* expected

```
uce-1732.nexus,609,4,94,609,41.37,193,724,358,570,591
uce-2349.nexus,576,0,123,576,35.36,0,508,287,324,609
uce-1985.nexus,601,0,121,601,37.16,0,549,345,325,584
uce-6550.nexus,531,0,147,531,37.1,0,501,274,317,501
uce-2539.nexus,551,0,79,551,40.17,0,378,331,333,611
uce-6187.nexus,540,0,97,540,41.54,0,492,333,340,455
uce-865.nexus,610,0,77,610,39.62,0,507,370,355,598
uce-2120.nexus,649,17,226,649,38.38,191,743,431,492,739
uce-4179.nexus,675,43,389,675,45.33,304,722,532,554,588
uce-1484.nexus,466,0,64,466,36.12,0,454,282,223,439
uce-7542.nexus,587,0,89,587,34.36,0,561,257,348,595
uce-6095.nexus,497,0,87,497,42.86,0,378,336,303,474
uce-3046.nexus,691,18,190,691,35.22,223,753,448,447,893
uce-553.nexus,782,3,161,782,41.56,277,841,601,584,825
uce-7014.nexus,490,17,175,490,42.94,211,465,402,349,533
```

## muscle-gblocks

```bash
python ../bin/align/phyluce_align_get_gblocks_trimmed_alignments_from_untrimmed \
    --alignments muscle \
    --output muscle-gblocks \
    --cores 4
```

* get stats:

```bash
python ../bin/align/phyluce_align_get_align_summary_data \
	--alignments muscle-gblocks \
	--output-stats muscle.gblocks.stats
```

* expected

```
uce-1732.nexus,609,4,93,609,41.44,194,723,358,571,590
uce-2349.nexus,652,0,125,652,35.62,78,557,328,341,652
uce-1985.nexus,729,0,139,729,37.97,122,619,392,392,662
uce-6550.nexus,711,0,166,711,37.17,180,635,361,365,592
uce-2539.nexus,629,0,84,629,40.8,78,411,369,369,660
uce-6187.nexus,591,0,101,591,41.52,51,519,358,357,488
uce-865.nexus,665,0,85,665,39.77,59,530,394,376,636
uce-2120.nexus,606,15,204,606,38.49,135,698,411,470,710
uce-4179.nexus,550,33,247,550,44.23,172,619,460,437,512
uce-1484.nexus,507,0,74,507,35.47,41,486,295,230,469
uce-7542.nexus,695,0,98,695,33.76,121,614,281,382,687
uce-6095.nexus,528,0,90,528,42.5,31,406,349,311,487
uce-3046.nexus,661,19,168,661,35.69,192,726,434,441,851
uce-553.nexus,781,3,152,781,41.16,284,859,588,581,812
uce-7014.nexus,464,13,143,464,41.89,185,446,375,325,525`
```

## muscle-trimal

```bash
python ../bin/align/phyluce_align_get_trimal_trimmed_alignments_from_untrimmed \
    --alignments muscle \
    --output muscle-trimal \
    --cores 4
```

```bash
python ../bin/align/phyluce_align_get_align_summary_data \
	--alignments muscle-trimal \
	--output-stats muscle.trimal.stats
```

* expected

```
aln,length,sites,differences,characters,gc content,gaps,a count, c count, g count, t count
uce-1732.nexus,609,4,94,609,41.37,193,724,358,570,591
uce-2349.nexus,574,0,113,574,35.48,0,512,292,319,599
uce-1985.nexus,609,0,116,609,37.6,0,548,344,343,592
uce-6550.nexus,540,0,141,540,36.85,0,511,286,311,512
uce-2539.nexus,551,0,75,551,40.35,0,379,332,335,607
uce-6187.nexus,540,0,92,540,41.42,0,488,325,346,461
uce-865.nexus,610,0,70,610,39.62,0,506,368,357,599
uce-2120.nexus,642,15,226,642,38.2,152,754,430,493,739
uce-4179.nexus,597,36,292,597,44.29,198,686,498,472,534
uce-1484.nexus,468,0,61,468,35.68,0,457,286,215,446
uce-7542.nexus,586,0,85,586,34.41,0,556,260,345,597
uce-6095.nexus,497,0,85,497,42.79,0,381,334,304,472
uce-3046.nexus,665,21,175,665,35.65,169,739,440,448,864
uce-553.nexus,794,3,166,794,41.15,289,866,601,587,833
uce-7014.nexus,487,15,162,487,41.65,200,472,388,340,548
```

# Remove locus Names

```bash
python ../bin/align/phyluce_align_remove_locus_name_from_nexus_lines \
    --alignments mafft-trimal \
    --output mafft-trimal-clean \
    --cores 4
```

```bash
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

* expected

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

* expected phylip header

```
4 9257
```