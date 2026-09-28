# MIMICRY: Quick Start Guide
## Introduction
 
Welcome to MIMICRY! 
MIMICRY is a web-based resource for exploring cancer-related molecular mimicry data. 
The database is a curation of neoantigens and tumour-associated antigens for MHC class I that share peptide sequences with microbial species of the human microbiome. An epitope with such a match is a candidate mimicry-epitope. 

This quick start guide introduces the main workflow for using MIMICRY and will help you move from a biological question to relevant database results.

- Access to MIMICRY: https://cgi.cse.unsw.edu.au/~z3021002/mimicry/
- Users: cancer researchers interested in immunogenic peptides  
 
## Search the Database
1. Open MIMICRY in your web browser.
2. From the landing page, you can search for a peptide, gene or microbe of interest within the database using the available text box in the filter sidebar.
3. The search will dynamically update the remaining pages.
4. Click on the "Peptide Search" tab to view all matching peptides.

### Tips for searching
- Start with a relatively broad query if you are exploring the database for the first time - ie pick a cancer type or microbe of interest.
- Add filters when the initial search returns too many results.
- Check identifiers carefully when searching for a particular gene, protein, peptide, cancer type, or other biological entity.
- Try removing optional filters if your query produces no results. 

## Inspect a Record 
1.  In the "Peptide Search" tab, you can look at the resulting data in more detail by clicking on the peptide.
2.  This option will allow users to evaluate an individual MIMICRY entry.
   
### Tips for interpreting the candidates 
- Identify the microbes involved in the mimicry relationship and their context and if they are relevant to your question. 
- Examine the binding score and HLA-alleles, and also the number of matches across different contexts.
- Check the associated cancer type, organism, experimental context, and supporting references where available.
_Important: A database match should be treated as a candidate for further investigation rather than, by itself, proof of a biological or clinical relationship._

## Download Data 
1. In the "Peptide Search" tab, you can download the resulting data table as a CSV file.
2. Select the columns you wish to use, and then click on "Download CSV'. 

## Further analysis
There are multiple additional tabs to work through more in-depth analyses. We link the individual guides below:  

[Cancer & Mutation](how-to-cancer.md) -> exploring the peptide data and its origins. 

[Microbiome and Body Locations](how-to-micro.md) -> highlighting the microbial mimicry landscape. 

[HLA Binding](how-to-bind.md) -> characterisation of the immunogenicity of the data. 

[Population Coverage](how-to-popcov.md) -> insights into the immune response  

[Custom Chart](how-to-custom.md)-> a sandbox for users  

### Need More Help?
If you are new to MIMICRY, start with a broad search and explore several individual records before applying more restrictive filters.
For publications or reports using MIMICRY, record the database version or access date and describe the search criteria used so that the analysis can be reproduced.
