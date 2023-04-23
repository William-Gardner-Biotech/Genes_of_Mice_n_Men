
<h1 style="color: Brown"><center>Genes of Mice and Men</center></h1>

![](https://i.imgur.com/ML1UCWo.png)

[TOC]

<h2 style="color: Brown"><center>Background</center></h1>


### Research Question:

<u>Which genes are conserved between <i>Homo sapiens</i> and <i>Mus musculus</i> and can be used to study human disease? Specifically, which genes can we study in a mouse model to research human cerebellum function and disease?</u>

Mice and humans share approximately 70% of the same protein-coding sequences (with protein coding sequences making up approximately 1.5% of each respective 3 billion base pair genome). This, along with fast reproduction, ease of care, and having short life spans, make mice a powerful test subject to use for furthering our understanding of complex human genetic diseases. Therefore, in order to maintain the progress made in research in gene annotation and novel drug discovery it is vital to find orthologous regions of DNA between mice and humans. Exploring differences between gene transcription across species may also help us to understand disease formation and mechanisms.[1][2]

Orthologous genes are the result of genes conserved between two species after a divergent event from a common ancestor. According to Nature[3] only 1% of human genes do not have a mouse ortholog, and between those orthologous genes is an average sequence similarity of 78.5%. If the goal is to utilize the mouse genome to make novel discoveries in realizing how a certain gene is playing its part in disease/function, it would be beneficial to know which genes share the highest sequence and tissue expression similarity. 

One of the anatomical areas that can be researched between humans and mice is the brain, more specifically the cerebellum, not only due to the aforementioned genetic similarities but also the anatomical resemblance. The cerebellum is located at the back of the head near where the spinal cord connects, and helps to coordinate and regulate a wide range of functions and processes in both your brain and body. [8] A disorder in the cerebellum can have effects such as; types of ataxia (failure of muscle control), brain degeneration, and cancer. The answer to what causes the development of such diseases/function may one day be found on the microscopic level of DNA, but first we must find the orthologs.

<h2 style="color: Brown"><center>Procedure</center></h1>

We took orthologous gene data from our own **BlastP** data table, Mouse and human brain transcription datatables from the **Human Protein Atlas** database. We then merged the three datasets by common gene name.

![](https://i.imgur.com/PhBwMry.png)


<h3 style="color: Orange"><center>Bash</center></h1>
The basic local alignment search tool is an algorithm designed to make sequence comparisons. We will use this tool to find significantly similar sequences within the human and mice proteomes.

1. Start by downloading the Mouse and Human proteomes from Uniprot. 
    (Link 10)
2. Using "mv" or "cp" command, move these downloads into a workable directory.
3. use gunzip to unzip the compressed fasta files
4. using blast+ program use the "makeblastdb" function and add in the unzipped human proteome
5. The mouse proteome is formatted the same as a multiFASTA file so we can use it as our query.
6. Using standard blast command line:
```
blastp -query mouse_proteome.fasta -db human_proteome.fasta \
    -out mouse_vs_human_blast_results.tab \
    -evalue 1e-20 \
    -outfmt 6\
    -max_target_seqs 1
```

Using -outfmt 6 we can get a much simpler data table that does not include all the additional comments like query number and number of hits.

Also using max_target_seqs 1 we reduce the amount of returned protein matches to only 1.

<h3 style="color: Orange"><center>Python</center></h1>

7. From this point we have a BLAST file but wanted only the gene names of the human gene and mouse gene. So using Python and the Regular Expression (re) module we were able to extract the names only and create our first data table.

![](https://i.imgur.com/8O52Erm.png)


```
import re
results_location = "/home/will/Proj_2/mouse_vs_human_blast_results.tab"

results = open(results_location, "r")
out_list = ["Mouse Gene, Human Gene"]
for line in results:
    splits = line.split("\t")
    mouse = (splits[0])
    subgroups_mouse = re.search("((\w*)$)", mouse)
    mouse = (subgroups_mouse.group(0))
    human = (splits[1])
    subgroups_human = re.search("((\w*)$)", human)
    human = (subgroups_human.group(0))
    addition = "\n"+mouse+","+human
    out_list.append(addition)

outfile = open("/home/will/Proj_2/Mouse_Human.csv", "w")

for i in out_list:
    outfile.write(i)

outfile.close(
```

![](https://i.imgur.com/7RMM29G.png)

*using the command: "blast_Mus_v_humo <- unique(Mus_v_humo)" we removed repition introduced by a problem created with the blastp argument, "max_target_seqs 1"

<h3 style="color: Orange"><center>RStudio</center></h1>

8. The second table we used was downloaded from ProteinAtlas. It contained human genes and tissue expression data expressed in nTPM (normalized Transcript Per Million) for humans. After removing a couple of unwanted data columns we were left with the following table:
![](https://i.imgur.com/Bg99xAk.png)

9. Using a column filter we removed all rows that did not include the "cerebellum" tissue. Now we had a table with nTPM readings only in the cerebellum and could merge it with our BLAST result table.
![](https://i.imgur.com/MS02BGk.png)

10. Using R, we merged the Blast results using "Human_Gene" and the human expression table "Gene_name." This returned the first merged table but to make meaningful comparisons between mice and human cerebllum genes we needed tissue level mouse RNA-seq data.

        merge_1 <- left_join(blast_Mus_v_humo, human_cereb, by=c("Human_gene"= "Gene_name"))
        
11. The table populated for all the human cerebellum genes and this is where our first level of tissue specificity is introduced. Many genes found using BLASTp were not expressed in the cerebellum so we removed their rows. Using two commands we removed genes that had no cerebellum (or very little <1) nTPM reads. 

        merge_1 %>% drop_na("Human_ENSG") 

        merge_1 <- filter(merge_1, Human_nTPM > 1)

12. Now using mouse data from Protein Atlas we filtered the table again for nTPM in the cerebllum and Mouse Gene name. The original data table was over 100 columns long but we used colnames to reduce this drastically. We included the Gene description and other general functional annotations for a more robust final data table.

![](https://i.imgur.com/OAkr3Zm.png)

13. Now using leftjoin we merged the mouse data table with our first merge table. We merged the mouse table onto the other table aligning on Mouse_gene and Gene. We then renamed the nTPM columns to be specific by species and moved them next to one another.

        merge_2 <- merge_2 %>% relocate(12, .after = 4)

14. Using the method from earlier we removed any N/A column that was a Mouse_gene with no associated annotation from the mouse table. What we noticed after this merge is that every mouse gene and human gene had the same gene name. This was very interesting because it showed us how genome annotations are created. Most genomes are referenced back to the human genome and named accordingly. Removing these N/A columns along with again filtering against genes with less than 1 nTPM returned our final table for species-species cerebellum specific gene expression table.

![](https://i.imgur.com/NFJ4R6S.png)

*Disease involvement is a column to the right.


```
setwd("C:/Users/Will/Desktop")

library(tidyverse)
library(dplyr)

Mus_v_humo <- read_csv(file = "Mouse_Human_clean.csv")

#view(Mus_v_humo)

mice <- read_tsv(file = "mouse_brain_category_rna_cerebral.tsv")

mouse_rough <- mice

#view(mouse_rough)

humans <- read_tsv(file = "rna_tissue_gtex.tsv")

humans = humans[-c(4,5)]

human_cereb <- filter(humans, Tissue == "cerebellum")

colnames(human_cereb)[4] = "Human_nTPM"

human_cereb <- human_cereb[-c(3)]

colnames(human_cereb)[1:2] = c("Human_ENSG", "Gene_name")

#view(human_cereb)

#view(human_cereb)

#view(Mus_v_humo)

blast_Mus_v_humo <- unique(Mus_v_humo)

colnames(blast_Mus_v_humo) = c("Mouse_gene", "Human_gene")

#view(blast_Mus_v_humo)

mice <- mice[c(1:17,93)]

mice_awk <- mice[-c(2,6,7,12:16)]

colnames(mice_awk)[5] = "Protein_desc"

view(mice_awk)

merge_1 <- left_join(blast_Mus_v_humo, human_cereb, by=c("Human_gene"= "Gene_name"))

merge_1 %>% drop_na("Human_ENSG") 

merge_1 <- filter(merge_1, Human_nTPM > 1)

view(merge_1)

merge_2 <- left_join(merge_1, mice_awk, by=c("Mouse_gene"= "Gene"))

merge_2 <- merge_2 %>% relocate(13, .after = 4)

merge_2 <- merge_2 %>% drop_na("Protein_desc") 

colnames(merge_2)[5] = "Mouse_nTPM"

colnames(merge_2)[13] = "RNA_tissue_specificity"

merge_2 <- merge_2[-c(3,6,8)]

view(merge_2)

write_csv( merge_2, file = "Final_table.csv")
```

An optional filtering to remove genes expressed throughout all tissues can be used to further refine the dataset to genes that may be more specific to the cerebllum. 

        merge_2 <- filter(merge_2, Mouse_nTPM > 1)

        merge_2 <- filter(merge_2, RNA_tissue_specificity != "Detected in all")
        
<h2 style="color: Brown"><center>Results/Data Utilization</center></h1>

![](https://i.imgur.com/t7cG6fN.png)

*Zoomed in Final data table.


The most applicable genes that could be used for research with mice would be the ones that show the highest sequence similarity and species to species expression. We have (1) identified orthologous proteins between mice and humans and (2) applied normalized transcription data for each gene from each respective species. Our data table holds over 3000 unique genes along with nTPM data between species and their respective annotations.

A researcher can look at our table and learn about the cerebellar function of genes such as: biological process, disease involvement, and molecular function and decide which genes to make perturbations to research a disease of interest. A quick analysis that can also be performed is to find which of the genes share a similar expression level to each other. Here is how we used our table to find genes whose transcription levels were within a user specified % similarity to each other using Python: 


```
results_location = "/[users director]/Final_table.csv"

results = open(results_location, "r")

header = results.readline()

# input your tolerance limit as a percentage difference

tolerance = input("Enter you percentage difference in similarity as a number:")

if tolerance.isnumeric():
    tolerance = float(tolerance)
    Similar_expression = []
    sp_sp_expression = open("/[user's directory]/Similar_spx2_" + str(tolerance) + ".csv", "w")
    sp_sp_expression.write(header)
    for line in results:
        splits = line.split(",")
        # [3] == human nTPM
        #print(splits[3], splits[4])
        human_nTPM = float(splits[3])
        mice_nPTM = float(splits[4])
        first_step = abs(human_nTPM-mice_nPTM)
        second_step = (mice_nPTM+human_nTPM)/2
        perc_diff = (first_step/second_step)*100
        if perc_diff < tolerance:
            sp_sp_expression.write(line)
        else:
            continue
    sp_sp_expression.close()
else:
    exit("Invalid input, Exiting program")
```

Another way that this table can be utilized is to make phylogenetic inferences on the evolution from the mouse cerebellum to the human cerebellum. Significant differences in expression of the same gene shows how the different species have diverged by illustrating how different proteins became more/less important for cerebellar processes for each respective species. We can also infer that genes with similar transcription levels between species are most likely more conserved in function.

<h2 style="color: Brown"><center>Discussion</center></h1>

We successfully merged transcriptional data across species to compare shared genes between humans and mice. We thought the utilization of BlastP to find orthologous genes would enable us to link together gene names that were not the same, however we found that all of the orthologous genes that had expression in cerebellum tissue already shared  a common gene name. We learned that the model used for creating gene names in novel genomes is heavily referenced back to the human genome. Which allows for easy comparisons. While the names of these genes were shared, our datatable shows that these genes were sometimes transcribed at significantly different levels between species. Our postulation is that the most fitting genes that can be used in studying human disease and function would be those that have similar transcription levels and gene interactions. Take for example the gene NOP56:

According to our table, NOP56 is a ribonucleoprotein that is involved with neurodegeneration and spinocerebellar ataxia. This gene showed similar cerebellum tissue expression level across species (Humans: 47nTPM, Mice: 39nTPM). In this way NOP56 could be utilized in research for these aforementioned diseases, but the interactions that this gene has with other genes will be needed to take into account to truly study the disease as a whole. If researchers wish to study NOP56 in humans using mice, they need to understand that the gene interactions and transcription levels could be different so they must adjust their research methods accordingly.

**Human NOP56 Gene Network:**
![](https://i.imgur.com/JrfSiVD.png)

**Mice NOP56 Gene Network:** 
![](https://i.imgur.com/g8UhAT3.png)

**Shared:**
* NOP58
* DKC1
* RRP9
* FBL
* NHP2L1
* NHP2

**Different:**
* MPHOSPH10 - human
* HEATR1 - human 
* KRR1 - human 
* WDR36 - human 
* WDR43 - mouse
* GM5616 - mouse
* FTSJ3 - mouse
* NIP7 - mouse

Moving into the future we would like to build this into a more robust database. A database that connects not just humans to mice but humans to other model organisms to guide disease research to delineate the best model organism. The database would be very reliant on RNA-Seq data and would be limited in scope due to the scarcity of this information source. We would also like to automate our process to create a database of similar tables that take into account user interest to produce a more specific/meaningful table. We can also expand the database past the cerebellum specific genes and have categories for every tissue or mutliple tissue genes. 

## References 

[1]Why Mouse matters
https://www.genome.gov/10001345/importance-of-mouse-genome#:~:text=On%20average%2C%20the%20protein%2Dcoding,they%20are%20required%20for%20function.

[2]Why use the mouse in research 
https://www.yourgenome.org/facts/why-use-the-mouse-in-research/

[3]Humanising the mouse genome
https://www.nature.com/articles/s41467-019-09716-7#:~:text=Approximately%201%25%20of%20human%20genes,numbers%20in%20the%20two%20species.

[4]Protein ID conversion tool - create name table 
https://www.genenames.org/download/custom/

[5]Mouse Proteome:
https://www.uniprot.org/proteomes/UP000000589

[6]Human Proteome:
https://www.uniprot.org/proteomes/UP000005640

[7]A comparison of human and mouse gene co-expression networks reveals conservation and divergence at the tissue, pathway and disease levels
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4654840/

[8]What is the cerebellum
https://my.clevelandclinic.org/health/body/23418-cerebellum

[9]Cerebellar Disorders
https://medlineplus.gov/cerebellardisorders.html

[10]Website Page containing our used databases
https://www.proteinatlas.org/about/download
Human [RNA GTEx brain region gene data]
Mice [RNA mouse brain region gene data]
