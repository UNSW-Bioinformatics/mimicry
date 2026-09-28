## How-to: population coverage analysis
The Population Coverage page answers for a given population, what share carries an HLA type that would present at least one of the mimicked peptides for a given cancer, tissue or gene. Population coverage is based on the distribution of HLA alleles in different populations. Because HLA allele frequencies vary geographically and between populations, a peptide predicted to interact with particular HLA alleles may have different estimated coverage across populations.

The calculation is the [IEDB method](https://tools.iedb.org/population/) and their allele frequencies, reimplemented in `bin/popcov.py`. 
Diploid genotypes are enumerated within each locus, loci are combined assuming independence, and where a locus's frequencies do not sum to 1
the shortfall becomes an untyped allele that binds nothing, which keeps the result a lower bound.
  
### Coverage calculator
From a MIMICRY result, open the Population Coverage tab to view the estimated population coverage associated with the peptide and its HLA predictions. One of the questions this can answer is whether the peptide is relevant across many populations, or if is its predicted HLA presentation concentrated in particular population. Higher coverage indicates that a greater proportion of individuals in that population are estimated to carry one or more of the relevant HLA alleles.

#### Population
Which population’s HLA frequencies to use. 
Coverage is entirely relative to this selection as an epitope set that reaches most of one population may reach far fewer of another, because HLA type varies enormously between them. Grouped by the geographic areas the IEDB data itself defines; only populations with usable class I frequencies are listed. 

#### Group by
Users can select to group by multiple fields, including: 
1. Cancer
2. Sample tissue
3. Gene
4. Antigen type
5. Variant type
6. Antigen source
7. Phylum
8. Genus
9. Body site

#### Epitopes per group
In addition, you can select how many peptides to consider per group. A vaccine carries a defined set, not every epitope that matched, and about 30 is the usual upper bound for an epitope-based T-cell vaccine. Epitopes are ranked by their best netMHCpan %rank, strongest first, with ties broken by how many alleles they bind. Coverage rises: the first few and then flattens.

#### Binding strength
As with filtering, binding strength allows the user to select peptides that are more likely to be true binders. 

#### Groups shown
This is a filter for the barplot displayed on the page, and in the tables below. 

### Outputs: 
1. A barplot of the groups selected and their fraction coverage for the selected population. This can be exported as an image or the data copied over to plot the same. 
2. A table of the peptides and what they bind is displayed below this. This can be used to explore the results in more detail. 
3. The coverage table displays the statistics around the coverage.
4. A final plot shows the peptides per individual in the population. This plot can also be exported or the data points copied over to plot the same. 
