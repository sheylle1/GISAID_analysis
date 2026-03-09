# Viral Sequence Analysis Pipeline for GISAID Submission

A Python pipeline for processing **Genome Detective outputs**, evaluating sequencing quality metrics, and preparing **FASTA files and reports for GISAID submission**.

The pipeline automatically parses `.assignments.json` files, evaluates genome coverage, generates submission-ready FASTA files, and produces summary reports for quality control.

---

# Features

- **Multi-virus support**
  - Influenza A
  - Influenza B
  - HIV
  - RSV
  - SARS-CoV-2  
  *(Easily extendable to additional viruses)*

- **Automatic parsing of Genome Detective outputs**
  - Supports both **old and new JSON formats**

- **Quality control for GISAID submission**
  - User-defined coverage thresholds
  - Depth of coverage filtering
  - Coverage percentage filtering

- **Influenza-specific analysis**
  - Segment-level analysis
  - HA (segment 4) and NA (segment 6) validation
  - Segment FASTA generation
  - Automatic segment concatenation

- **FASTA file generation**
  - Customizable FASTA headers
  - Optional multi-lab formatting
  - Automatic extraction of consensus sequences

- **Lab-aware metadata**
  - Supports **single or multiple labs**
  - Interactive assignment of samples to labs

- **Robust sample handling**
  - Automatic extraction of unique IDs (`C0XXXX` or `K0XXXX`)
  - Prevents duplicate processing
  - Handles control samples separately

- **Automated reporting**
  - GISAID eligibility summaries
  - Segment quality metrics
  - Full genome statistics

---

# Supported Viruses

| Virus | Genome Type | Special Processing |
|------|------|------|
| Influenza A | Segmented | HA/NA submission checks |
| Influenza B | Segmented | HA/NA submission checks |
| HIV | Non-segmented | Whole genome evaluation |
| RSV | Non-segmented | Whole genome evaluation |
| SARS-CoV-2 | Non-segmented | Whole genome evaluation |

Additional viruses can be added easily by modifying the `VIRUS_CONFIG` dictionary.

---

# Installation

```bash
git clone https://github.com/yourusername/viral-gisaid-pipeline.git
cd viral-gisaid-pipeline

pip install -r requirements.txt
```

---

# Usage

Basic command:

```bash
python Gisaid_analysis.py /path/to/gd_processed/
```

Example:

```bash
python Gisaid_analysis.py /analyses/CERI/.../gd_processed/
```

The input folder should contain **Genome Detective processed sample folders**.

---

# Pipeline Workflow

1. Copy `.assignments.json` files into a temporary directory.
2. Detect the selected virus in the samples.
3. Evaluate sequencing quality:
   - Depth of coverage
   - Genome coverage percentage
4. Determine **GISAID submission eligibility**.
5. Generate **summary reports**.
6. Extract consensus sequences from alignment files.
7. Generate **submission-ready FASTA files**.
8. Concatenate segments/genomes.
9. Produce combined FASTA datasets for downstream analysis.

---

# Input Structure

The pipeline expects a directory structured as follows:

```
gd_processed/
├── Sample_1/
│   ├── sample.assignments.json
│   ├── *.alignment-nt.fasta
│
├── Sample_2/
│   ├── sample.assignments.json
│   ├── *.alignment-nt.fasta
```

Each sample folder should contain:

- `.assignments.json` file from **Genome Detective**
- Alignment FASTA files (`*alignment-nt.fasta`)

---

# Output Structure

```
project/

├── GD_Assignments/        # Temporary JSON files

├── Outputs/               # Generated FASTA files
│   └── {LimsID}/
│       ├── segment1.fasta
│       ├── segment2.fasta
│       └── {LimsID}_all_segments.fasta

├── Summary_files/         # Quality reports

└── viral_analysis.log     # Processing log
```

---

# Generated Reports

| File | Description |
|-----|-----|
| `gisaid_submission_status_{virus}.csv` | Sample-level submission eligibility |
| `influenza_{type}_segments.csv` | Segment-level metrics (Influenza only) |
| `full_genome_info_{virus}.csv` | Whole genome sequencing statistics |

---

# FASTA Outputs

### Influenza
- Individual segment FASTA files
- Concatenated segment FASTA per sample
- Combined FASTA file for all samples

### Non-segmented viruses
- Single FASTA per sample
- Combined FASTA dataset

Headers can be customized using placeholders:

```
<LimsID>
<gene>
<lab>
<isolate>
```

Example header:

```
>A/South Africa/NHLS-CERI-C01234/2026_HA
```

---

# Quality Thresholds

User-defined thresholds determine **GISAID eligibility**.

Default values:

| Metric | Default |
|------|------|
| Minimum depth | 10× |
| Minimum genome coverage | 80% |

These can be modified during runtime.

---

# Lab Assignment

The pipeline supports:

### Single-lab datasets
Assign all samples to one lab.

### Multi-lab datasets

Samples can be assigned interactively:

```
UCT K001-K050
NHLS C0100-C0200
```

---

# Logging

Processing steps are recorded in:

```
viral_analysis.log
```

Logs include:

- Sample processing
- File creation
- Missing alignment files
- Quality evaluation results
- Errors and warnings

---

# Extending the Pipeline

New viruses can be added by editing:

```python
VIRUS_CONFIG
```

Example structure:

```python
'virus_name': {
    'name': 'Virus Name',
    'segmented': False,
    'taxonomy_keywords': ['taxonomy name'],
    'common_name': 'VirusShortName'
}
```

---

# Requirements

- Python **3.9+**

Python packages:

- pandas
- pathlib
- json
- logging

Install with:

```bash
pip install -r requirements.txt
```

---

# Contributing

1. Fork the repository

2. Create a feature branch

```bash
git checkout -b feature/new-feature
```

3. Commit your changes

```bash
git commit -m "Add new feature"
```

4. Push to GitHub

```bash
git push origin feature/new-feature
```

5. Open a Pull Request

---
