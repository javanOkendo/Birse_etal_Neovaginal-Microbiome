# Birse_etal_Neovaginal-Microbiome
The neovaginal microbiome of transgender women post-gender reassignment surgery

January 17, 2020

These scripts detail the data processing steps utilized for the analysis of the 16S rRNA gene amplicon sequence data and metaproteome data obtained from transgender women's neovaginal and rectal samples as well as from Canadian and Swedish cisgender women's vaginal samples utilized in this study.

###Inventory:

##Metaproteome:

Datasets: 
"Bacterial proteins.csv"
"Human proteins_NV vs cis V_progenesis output.csv"
Codebook:
"Metaproteome R analysis codebook for data processing.Rmd"
Functions: 
"bacteria_pplot_clrs_v2.R"
"PPFunctions_noCLRS_edited theme_KB_shanhedit.R"

##16S rRNA gene analysis:

#Transgender Women:

Datasets:
"Neotrans_phyloseq_object.rds"
Metadata:
"Neotrans_Mapping.011620.txt"
"raw_abundance.phyloseq_study ids.txt"
Codebooks/Scripts:
"DADA2_101117.R"
"run_phyloseq_Neovag.R"
Functions:
"utils.R"

#Cisgender Women:

Datasets:
"Canadian cis vaginal 16S.raw.mothur.output.csv"
"Swedish 16S cis vaginal otu table.csv"
"Swedish cis vaginal V15L otu table.csv"
Codebook:
"Canadian cis vaginal_16S analysis_R Markdown_Jan.14.2020.Rmd"
"batch_file_all_commands_Microbiome manuscript analysis_01.22_final.sh"
Functions:
"Neovagina Paper_16S analysis_functions_Jan.14.2020"

