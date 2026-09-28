## How-to: population coverage analysis
The Population Coverage page answers for a given population, what share carries an HLA type that would present at least one of the mimicked peptides for a given cancer, tissue or gene.

The calculation is the [IEDB method](), reimplemented in `site/popcov.py` using only the standard library. 
Diploid genotypes are enumerated within each locus, loci are combined assuming independence, and where a locus's frequencies do not sum to 1
the shortfall becomes an untyped allele that binds nothing, which keeps the result a lower bound.
  
 
