## Simple code to run SCENIC on your single cell data

The fundamental big plus of single cell seq experiments is exploring the heterogeneity of a biological sample. 
This is particularly interesting (I think) if you are focusing on a single tissue type as it undergoes change
This might be cells in a developing tissue, in culture or cells in a disease process (maybe compared to normal tissue)

When you are doing the above, scRNAsq is essentially showing you the various _transcriptional states_ of a tissue. 
How many states exist within your sample is a tricky question, and kind of centres around some sort of combination
of how diverse the sample actually is, combined with your clustering resolution and how confident you are that each
cluster you annotate as a seperate state/type is actually that. 

Once you've identified these changing cell states, it makes sense to ask what is driving or regulating the changes.
In other words, **Which transcription factors are pulling the strings?**

TO keep it SCimple, SCENIC does two things:

**1)** It finds out which genes are co-expressed with TFs in the data, and writes that down in a table
**2)** It scrutinises the results to try and reduce noise (chance co-expression)

The latter step is more complex, but basically it uses databases (you need to download these and point SCENIC in their direction)
to check if the appropriate **transcription factor binding motifs** are significantly over-represented near the **transcription start site (TSS)**
of the genes picked out in step 1). This leaves you with a set of TFs and their targets, known together as a **'regulon'**



```markdown
Syntax highlighted code block

# Header 1
## Header 2
### Header 3

- Bulleted
- List

1. Numbered
2. List

**Bold** and _Italic_ and `Code` text

[Link](url) and ![Image](src)
```

For more details see [GitHub Flavored Markdown](https://guides.github.com/features/mastering-markdown/).

### Jekyll Themes

Your Pages site will use the layout and styles from the Jekyll theme you have selected in your [repository settings](https://github.com/SCimpleton/simpleSCENIC/settings). The name of this theme is saved in the Jekyll `_config.yml` configuration file.

### Support or Contact

Having trouble with Pages? Check out our [documentation](https://docs.github.com/categories/github-pages-basics/) or [contact support](https://support.github.com/contact) and we’ll help you sort it out.
