# OSC Genome UCSC Assembly Hub Pipeline

Pipeline for creating a UCSC Genome Browser Assembly Hub from the OSC genome, including annotation tracks, repeat tracks, and other genomic features.

Part of the **Handler et al., 2025** publication:

**The Drosophila OSC Genome: A Resource for Studies of Transposon and piRNA Biology**

## Overview

This repository contains the complete workflow for generating a UCSC Genome Browser Assembly Hub for the OSC genome. The pipeline creates a comprehensive browser hub with multiple annotation tracks including gene models, transposable elements, piRNA clusters, and other genomic features. The resulting hub enables interactive visualization and exploration of the OSC genome and associated datasets.

**Access the OSC Genome Browser Hub:**
https://genome-euro.ucsc.edu/s/Brennecke%2DLab/OSC_r1.01_Handler_et.al._2025

## Repository Structure

```
├── script-files/           # Core hub generation scripts
├── utility-files/          # Helper scripts and configuration files
├── .gitignore             # Git ignore file
└── create_assmblyHub.sh   # Main submission script for hub creation
```

## Pipeline Components

### Genome Preparation
Scripts for formatting the genome assembly for UCSC Browser compatibility, including sequence indexing and chromosome naming.

### Track Generation
Tools for creating and formatting various annotation tracks:
- Gene annotations
- Transposable element annotations
- piRNA cluster annotations
- Repeat masking tracks
- Custom data tracks

### Hub Configuration
Scripts for generating hub configuration files (`hub.txt`, `genomes.txt`, `trackDb.txt`) and organizing track hierarchies.

### Data Formatting
Conversion of annotation files to UCSC-compatible formats (BED, bigBed, bigWig, BAM, etc.).

## Requirements

- Apptainer
- Flybase reference files placed into the utility-files directory (required as automatic downloads do not work properly anymore)
  -   dmel-all-CDS-r6.*.fasta.gz
  -   dmel-all-miRNA-r6.*.fasta.gz
  -   dmel-all-miscRNA-r6.*.fasta.gz
  -   dmel-all-ncRNA-r6.*.fasta.gz
  -   dmel-all-pseudogene-r6.*.fasta.gz
  -   dmel-all-r6.*.gff.gz
  -   dmel-all-transcript-r6.*.fasta.gz
  -   dmel-all-tRNA-r6.*.fasta.gz
  -   fbgn_annotation_ID_fb_*.tsv.gz


## Usage

Run the main hub creation pipeline using:

```bash
bash create_assmblyHub.sh
```

Modify configuration parameters in the script files to customize track settings, colors, and display options.

## Output

The pipeline produces:
- Complete UCSC Assembly Hub directory structure
- Indexed genome sequence (`.2bit` format)
- Formatted annotation tracks (bigBed, bigWig, BAM, etc.)
- Hub configuration files
- Track database files


## Related Resources

### Main Publication Repository
https://github.com/BrenneckeLab/Handler_2025-OSC-genome


## Citation

Please find the proper citation in https://github.com/BrenneckeLab/Handler_2025-OSC-genome

## Contact

For questions or additional information, please contact:
dominik.handler@imba.oeaw.ac.at

## License

MIT License
