---
title: Population coverage in MIMICRY
shortTitle: <subject> # Max 31 characters
intro: 'Article intro. See tips for a great intro below.'
product: "{{ optional product callout }}"
contentType: how-tos
versions:
  - 1 
---

The Population Coverage page answers: for a given population, what share carries
an HLA type that would present at least one of the mimicked peptides for a given
cancer, tissue or gene.

The calculation is the IEDB method, reimplemented in `site/popcov.py` using only
the standard library. The IEDB distribution needs numpy and ships its
frequencies as a pickle, neither of which survives a CGI process on CSE. The
maths is unchanged: diploid genotypes are enumerated within each locus, loci are
combined assuming independence, and where a locus's frequencies do not sum to 1
the shortfall becomes an untyped allele that binds nothing, which keeps the
result a lower bound.
 


{% comment %}
Follow the guidelines in https://docs.github.com/contributing/writing-for-github-docs/content-model to write this article.
Great intros give readers a quick understanding of what's in the article, so they can tell whether it's relevant to them before moving ahead. For more tips, see https://docs.github.com/contributing/writing-for-github-docs/content-model
For product callout info, see https://github.com/github/docs/tree/main/content#product
For product version instructions, see https://github.com/github/docs/tree/main/content#versioning
Remove these comments from your article file when you're done writing
{% endcomment %}

## Procedural section header here

{% comment %}
Include prerequisite information or specific permissions information here.
Then write procedural steps following the instructions in https://docs.github.com/contributing/style-guide-and-content-model/style-guide#procedural-steps.
Check if there's already a reusable string for the step you want to write in https://github.com/github/docs/tree/main/data/reusables. Look at the source file for a procedure located in the same area of the user interface to find reusables.
{% endcomment %}

## Optionally, another procedural section here

{% comment %}
Keep adding procedures until you've finished writing your article.
{% endcomment %}

## Further reading

{% comment %}
Optionally, include a bulleted list of related articles the user can reference to extend the concepts covered in this article. Consider linking to procedural articles or tutorials that help the user use the information in your article.
{% endcomment %}

- [Article title](article-URL)
