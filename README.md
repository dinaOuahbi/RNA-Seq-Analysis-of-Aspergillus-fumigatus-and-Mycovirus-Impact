# RNA-Seq Analysis of *Aspergillus fumigatus* and Mycovirus Impact

## Overview
This repository provides a workflow for processing and analyzing RNA-Seq data of *Aspergillus fumigatus* with a focus on differential expression analysis. The goal is to investigate the presence of mycoviruses and their impact on gene expression, pathogenicity, and antifungal resistance.

## Background
### *Aspergillus fumigatus* and Its Clinical Relevance
*Aspergillus fumigatus* is a filamentous fungus responsible for severe respiratory diseases, particularly in immunocompromised individuals. It is a major causative agent of invasive aspergillosis (IA), a life-threatening infection with high mortality rates. The pathogen exhibits resistance to antifungal treatments, making it crucial to explore novel factors influencing its virulence.

### Mycoviruses and Their Role
Mycoviruses are RNA or DNA viruses that infect fungi without causing apparent cytopathic effects. Recent studies suggest that mycoviruses can modulate fungal virulence, stress response, and secondary metabolism. In *A. fumigatus*, certain mycoviruses have been associated with reduced growth, altered pathogenicity, and potential attenuation of virulence, making them intriguing targets for fungal disease control.

## Objectives
- Process and analyze RNA-Seq data from *A. fumigatus* strains with and without mycoviral infection.
- Perform differential gene expression analysis to identify key pathways affected by mycovirus presence.
- Investigate gene expression patterns related to virulence, antifungal resistance, and immune interactions.
- Generate visualizations and statistical summaries of the findings.

## Pipeline Overview
1. **Raw Data Preprocessing**
   - Quality control using FastQC
   - Adapter trimming with Trimmomatic
   
2. **Read Alignment**
   - Mapping to the *A. fumigatus* reference genome using HISAT2 or STAR
   
3. **Quantification**
   - Gene expression quantification using featureCounts or Salmon
   
4. **Differential Expression Analysis**
   - Performed with DESeq2 or edgeR
   - Identification of differentially expressed genes (DEGs) between infected and non-infected strains
   
5. **Functional Annotation & Pathway Analysis**
   - GO and KEGG enrichment analysis
   - Identification of pathways related to virulence and stress response
   
6. **Data Visualization**
   - Volcano plots, heatmaps, and PCA plots for comparative analysis
   
## Requirements
- Python (>=3.8)
- R (>=4.0)
- FastQC
- Trimmomatic
- HISAT2 / STAR
- featureCounts / Salmon
- DESeq2 / edgeR
- ggplot2, pheatmap (for visualization)

## Usage
### 1. Clone the repository
```bash
git clone [https://github.com/dinaOuahbi/mycovirus.git]
cd mycovirus
```

### 2. Run the pipeline
Ensure dependencies are installed and raw FASTQ files are in the `data/` directory.
```bash
sbatsh launch.sh
```
## Results & Interpretation
- **Differentially expressed genes**: Identifies key genes influenced by mycovirus infection.
- **Pathway enrichment**: Highlights functional pathways affected in *A. fumigatus*.
- **Impact on virulence**: Insights into how mycovirus presence modulates pathogenic potential.

## in
> [https://www.linkedin.com/in/dina-ouahbi-963a56338/]




