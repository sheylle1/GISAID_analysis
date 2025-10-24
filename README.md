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

## 🧾 Outputs

The pipeline generates the following structure and files:

### 📂 Directory Structure
```text
project/
├── GD_Assignments/          # Temporary JSON files
├── Outputs/                 # Final FASTA files
│   └── {LimsID}/            # Per-sample folders
├── Summary_files/           # CSV reports
└── viral_analysis.log       # Processing log
```

## Report Files

| File | Description |
|------|-------------|
| `gisaid_submission_status_{virus}.csv` | GISAID eligibility summary |
| `influenza_{type}_segments.csv` | Detailed segment information (Influenza only) |
| `full_genome_info_{virus}.csv` | Complete genome metrics |

## FASTA Files

- Individual segment files (for influenza)
- Concatenated genome files  
- Custom FASTA headers include lab and sample information

## Configuration

**Coverage Thresholds**
- Minimum depth of coverage: 10x (default)
- Minimum genome coverage percentage: 80% (default)

**Virus Configuration**
- Easily extendable in the `VIRUS_CONFIG` dictionary to support new viruses

## Requirements

- Python 3.9+
- See `requirements.txt` for full dependency list

## Contributing

1. Fork the repository
2. Create your feature branch:
   ```bash
   git checkout -b feature/amazing-feature
   ```
3. Commit your changes:
   ```bash
   git commit -m "Add amazing feature"
   ```
4. Push to your branch:
   ```bash
   git push origin feature/amazing-feature
   ```
6. 
