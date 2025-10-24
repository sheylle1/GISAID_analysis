# Viral Sequence Analysis Pipeline for GISAID Submission

A Python pipeline for processing Genome Detective outputs, evaluating quality metrics, and preparing files for GISAID submission.

## Features

- **Multi-virus support**: Influenza (A/B), HIV, RSV, and easily extensible to others
- **Quality evaluation**: User-defined coverage thresholds for GISAID eligibility
- **Automated FASTA generation**: Customizable headers with lab information
- **Segmented virus handling**: Specialized processing for influenza segments
- **Batch processing**: Handles multiple samples efficiently
- **Comprehensive reporting**: CSV summaries and quality control metrics

## Supported Viruses

- **Influenza A & B**: Segmented genome analysis with HA/NA quality checks
- **HIV**: Single-segment analysis
- **RSV**: Single-segment analysis
- **Extensible**: Easy to add new viruses via configuration

## Installation

```bash
git clone https://github.com/yourusername/viral-gisaid-pipeline.git
cd viral-gisaid-pipeline
pip install -r requirements.txt
```

## Usage

### Basic Usage
```bash
python Gisaid_analysis.py /path/to/gd_processed/folders
```
## Input Structure

The script expects a folder containing subdirectories for each sample, with:

.assignments.json files (Genome Detective outputs)

Alignment FASTA files

## Outputs 

The pipeline generates:

Directory Structure
text
project/
├── GD_Assignments/          # Temporary JSON files
├── Outputs/                 # Final FASTA files
│   └── {LimsID}/           # Per-sample folders
├── Summary_files/           # CSV reports
└── viral_analysis.log      # Processing log

## Report Files
gisaid_submission_status_{virus}.csv: GISAID eligibility summary

influenza_{type}_segments.csv: Detailed segment information (influenza only)

full_genome_info_{virus}.csv: Complete genome metrics

## FASTA Files
Individual segment files (influenza)

Concatenated genome files

Custom headers with lab/sample information

## Configuration

### Coverage Thresholds
- Minimum depth of coverage (default: 10x)
- Minimum genome coverage percentage (default: 80%)

### Virus Configuration
Easily extendable in VIRUS_CONFIG dictionary for new viruses.

## Requirements
Python 3.9+
See requirements.txt for dependencies

## Contributing

1. Fork the repository
2. Create a feature branch (git checkout -b feature/amazing-feature)
3. Commit your changes (git commit -m 'Add amazing feature')
4. Push to the branch (git push origin feature/amazing-feature)
4. Open a Pull Request