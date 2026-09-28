## How to: Explore HLA binding
The HLA binding tab lets you explore predicted binding between MIMICRY peptides and HLA class I alleles. The user can use this tab to identify peptides with stronger predicted HLA binding, determine which HLA alleles bind the largest numbers of peptides, and compare binding predictions across peptide–HLA combinations.

The tab contains four main components:
### Setting the binding thresholds
1. At the top of the tab, Where the binding threshold sits allows you to define the %rank cut-offs used to classify peptides as strong or weak binders. Lower %rank values indicate stronger predicted binding.
2. The default thresholds are:
- Strong binder: %rank ≤ 0.5
- Weak binder: %rank ≤ 2
3. Try different thresholds if you want to explore how more or less stringent definitions affect the results.

### Binding-rank distribution
The Binding-rank distribution shows the distribution of the best predicted HLA-binding score for each peptide.
Peptides towards the left-hand side have lower %rank values and therefore stronger predicted HLA binding.

### HLA allele coverage
The HLA allele coverage chart compares the number of peptides associated with different HLA alleles.
Each horizontal bar represents an HLA allele. The displayed bar distinguishes peptides meeting the strong and weak binding criteria.
Longer bars indicate that more peptides meet the relevant binding criterion for that HLA allele.
Strong binders are overlayed. 

### Peptide × HLA binding heatmap
The Peptide × HLA binding heatmap provides a detailed view of HLA-binding predictions, with the top 20 displayed.
Rows represent individual peptide sequences.
Columns represent HLA alleles.
Each coloured cell represents the predicted binding result for a particular peptide–HLA pair.
The colour scale corresponds to the displayed binding %rank. 
