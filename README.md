[![Tests](https://github.com/faircloth-lab/phyluce/actions/workflows/main.yml/badge.svg)](https://github.com/faircloth-lab/phyluce/actions/workflows/main.yml)

# phyluce: software for UCE (and general) phylogenomics

[phyluce][1] (phy-**loo**-chee) is a software package that was initially developed for analyzing data collected from ultraconserved elements in organismal genomes.

The package now includes a number of tools spanning:

* the assembly of raw read data to contigs
* the identification of UCE contigs
* parallel alignment generation, alignment trimming, and alignment data summary
  methods in preparation for analysis
* alignment and SNP calling using UCE or other types of raw-read data.

As it stands, the [phyluce][1] package is useful for analyzing both data collected from UCE loci and also data collection from other types of loci for phylogenomic studies at the species, population, and individual levels.

Please see the [Documentation][2] for additional information on installing and using [phyluce][1]

License
-------

3-clause BSD. See [License][3] for more information.

Contributions
--------------

[phyluce][1] is open-source (see [License][3]), and we welcome contributions from anyone who is interested in contributing.  To contribute code, please make a pull request on github_.  The issue tracker for phyluce_ is also available [on github][4].

Issues
------

If you have an issue, please ensure that you are experiencing this issue on a supported OS (see [Installation][5]) using the [conda][6] installation of [phyluce][1].  If possible, **please submit a test case demonstrating the issue** and indicate which platform and which phyluce version you are using.

Citing
------

If you use the [phyluce][1] code in any form, please cite the following manuscripts:

> Faircloth BC. 2015. PHYLUCE is a software package for the analysis of conserved genomic loci.  Bioinformatics. doi: [10.1093/bioinformatics/btv646](https://doi.org/10.1093/bioinformatics/btv646).

If you are processing UCE data that you have collected by targeted enrichment using our probes/protocols, please also cite the following manuscript, which describes the first use of the general approach:

> BC Faircloth, McCormack JE, Crawford NG, Harvey MG, Brumfield RT, Glenn TC. 2012. Ultraconserved elements anchor thousands of genetic markers spanning multiple evolutionary timescales. Systematic Biology 61: 717–726. doi: [10.1093/sysbio/SYS004](http://doi.org/10.1093/sysbio/SYS004).

[1]: https://github.com/faircloth-lab/phyluce "Link to this repository"
[2]: http://phyluce.readthedocs.io/ "Link to the documentation"
[3]: https://github.com/faircloth-lab/phyluce/blob/main/LICENSE "Link to the LICENSE"
[4]: https://github.com/faircloth-lab/phyluce/issues "Link to phyluce ISSUES"
[5]: http://phyluce.readthedocs.org/en/installation.html "Link to Installation"
[6]: https://docs.conda.io/en/latest/ "Link to conda documentation"