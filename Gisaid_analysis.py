#!/usr/bin/env python3.9
"""
Viral Sequence Analysis Pipeline for GISAID Submission

This script processes GD processed sequencing 
data evaluates quality metrics, and prepares files 
for GISAID submission.

Usage:
python /path/to/Gisaid_analysis.py /path/to/gd_processed/folders
Example:
python Gisaid_analysis.py /analyses/CERI/CERI0025.RPS004r.FluREP.VCS.b06P2_100625/gd_processed/

The script expects the path to the folder containing all GD-processed sample folders.

Features:
- Quality evaluation based on user-defined coverage thresholds
- Automated FASTA file generation with customizable headers
- Summary reports for quality control
- Support for multiple labs/sources
"""

import os
import json
import logging
from pathlib import Path
import pandas as pd
import re
import glob

# Virus configuration - makes it easy to add new viruses
VIRUS_CONFIG = {
    'influenza': {
        'name': 'Influenza',
        'segmented': True,
        'types': {
            'A': {
                'taxonomy_keywords': ['alphainfluenzavirus influenzae'],
                'segments': {
                    1: 'PB2', 2: 'PB1', 3: 'PA', 4: 'HA', 
                    5: 'NP', 6: 'NA', 7: 'MP', 8: 'NS'
                }
            },
            'B': {
                'taxonomy_keywords': ['betainfluenzavirus influenzae'],
                'segments': {
                    1: 'PB1', 2: 'PB2', 3: 'PA', 4: 'HA', 
                    5: 'NP', 6: 'NA', 7: 'MP', 8: 'NS'
                }
            }
        }
    },
    'hiv': {
        'name': 'HIV',
        'segmented': False,
        'taxonomy_keywords': ['human immunodeficiency virus'],
        'common_name': 'HIV'
    },
    'rsv': {
        'name': 'RSV',
        'segmented': False,
        'taxonomy_keywords': ['orthopneumovirus hominis'],
        'common_name': 'RSV'
    }
}

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('viral_analysis.log'),
        logging.StreamHandler()
    ]
)

def normalize_limsid(raw_id):
    """
    Normalize a raw LIMS ID by removing suffixes for non-control samples.
    
    Parameters:
        raw_id (str): The raw sample ID from the filename.

    Returns:
        str: Normalized LIMS ID.
    """
    if "PC" in raw_id or "NC" in raw_id or "Neg" in raw_id or "Pos" in raw_id:
        return raw_id  # Keep controls as-is
    elif "-R" in raw_id or "-r" in raw_id:
        return raw_id #These are repeats
    elif "_" in raw_id :
        return raw_id.split('_')[0]
    else:
        return raw_id.split('-')[0]

def create_temp_folder(input_folder):
    temp_folder = os.path.join(os.getcwd(), 'GD_Assignments')
    os.makedirs(temp_folder, exist_ok=True)
    logging.info(f"Created temp folder: {temp_folder}")
    
    for root, _, files in os.walk(input_folder):
        for file in files:
            if file.endswith('.assignments.json'):
                src_path = os.path.join(root, file)
                dest_path = os.path.join(temp_folder, file)
                if os.path.abspath(src_path) != os.path.abspath(dest_path):
                    try:
                        with open(src_path, 'rb') as src_file:
                            data = src_file.read()
                            if not data:
                                continue
                            else:
                                with open(dest_path, 'wb') as dest_file:
                                    dest_file.write(data)
                    except Exception as e:
                        logging.error(f"Failed to copy {src_path}: {e}")
    return temp_folder

def get_virus_type():
    """
    Prompt the user to select the virus type from a predefined list (e.g., Influenza, HIV, RSV).
    
    For Influenza, a follow-up prompt asks for subtype A/B/C.

    Returns:
        tuple: (virus_type (str), virus_key (str))
            - virus_type: User-selected virus type or subtype (e.g., 'A', 'HIV')
            - virus_key: Key for looking up metadata in the VIRUS_CONFIG dictionary
    """
    # Display menu of available viruses
    print("\nSelect virus type:")
    for i, (virus_key, virus_data) in enumerate(VIRUS_CONFIG.items(), 1):
        print(f"{i}. {virus_data['name']}")
    
    # Get and validate user input
    choice = input(f"Enter your choice (1-{len(VIRUS_CONFIG)}): ").strip()
    while not choice.isdigit() or int(choice) not in range(1, len(VIRUS_CONFIG) + 1):
        print("Invalid choice. Please try again.")
        choice = input(f"Enter your choice (1-{len(VIRUS_CONFIG)}): ").strip()
    
    # Resolve virus key from numeric selection
    virus_key = list(VIRUS_CONFIG.keys())[int(choice) - 1]
    
    # If Influenza, follow up with type-specific prompt (A/B/C)
    if virus_key == 'influenza':
        return get_influenza_type(), virus_key
    
    return virus_key.upper(), virus_key

def get_influenza_type():
    """
    Prompt the user to select which type of influenza virus samples are being analyzed.

    Options:
        A - Influenza A
        B - Influenza B
        C - Both A and B (combined analysis)

    Returns:
        str: One of 'A', 'B', or 'C' representing the influenza type selected.
    """
    print("\nWhat type of influenza samples are you looking at?")
    print("A. Influenza A")
    print("B. Influenza B")
    print("C. Both")
    
    # Prompt user input and validate
    choice = input("Enter your choice (A/B/C): ").upper().strip()
    while choice not in ['A', 'B', 'C']:
        print("Invalid choice. Please enter A, B, or C.")
        choice = input("Enter your choice (A/B/C): ").upper().strip()
    
    return choice

def get_coverage_thresholds():
    """
    Prompt the user to define quality thresholds for GISAID submission.

    Users may define:
        - Minimum depth of coverage (e.g., 10x)
        - Minimum genome coverage percentage (e.g., 80%)

    If left blank or invalid, defaults are applied (depth=10.0, coverage=80.0).

    Returns:
        tuple:
            float: minimum depth of coverage
            float: minimum genome coverage percentage
    """
    print("\nSet coverage thresholds for GISAID submission:")
    min_depth = input("Minimum depth of coverage (default 10): ").strip()
    min_cov = input("Minimum coverage percentage (default 80): ").strip()
    
    try:
        # Convert input to float or use defaults if blank
        min_depth = float(min_depth) if min_depth else 10.0
        min_cov = float(min_cov) if min_cov else 80.0
    except ValueError:
        # Handle non-numeric inputs gracefully
        logging.warning("Invalid input. Using default thresholds (depth=10, cov=80%)")
        min_depth, min_cov = 10.0, 80.0
    
    return min_depth, min_cov

def parse_json_file(json_file, virus_type, virus_key, min_depth, min_cov):
    """
    Parse a Genome Detective JSON `.assignments.json` file and extract viral information 
    relevant for downstream analysis and GISAID submission.

    Handles both:
        - Segmented viruses (e.g., Influenza A/B)
        - Non-segmented viruses (e.g., HIV, RSV)

    Parameters:
        json_file (str): Path to the input JSON file.
        virus_type (str): Specific virus subtype (e.g., 'A', 'B', or 'C' for influenza).
        virus_key (str): Key used to look up virus metadata in the VIRUS_CONFIG dictionary.
        min_depth (float): Minimum depth of coverage threshold for inclusion.
        min_cov (float): Minimum genome coverage (%) threshold for inclusion.

    Returns:
        list[dict] | None: A list of result dictionaries per valid strain, or None if no matches found.
    """

    try:
        # First check if file is empty
        if os.path.getsize(json_file) == 0:
            logging.warning(f"Empty file detected: {json_file}")
            return None
            
        with open(json_file, 'r') as f:
            try:
                data = json.load(f)
            except json.JSONDecodeError as e:
                logging.error(f"Invalid JSON in {json_file}: {e}")
                return None
            

        # Normalize the LIMS ID (remove suffixes unless it's a control)
        lims_id = normalize_limsid(os.path.basename(json_file).split('.')[0])

        # Get list of all strains from the file
        strains = data.get('data', {}).get('attributes', {}).get('strains', [])

        results = []

        for strain in strains:
            taxonomy = strain.get('taxonomyName', '').lower()

            # --- Match strain to correct virus type ---
            if virus_key == 'influenza':
                # Influenza-specific matching based on subtype
                if (
                    (virus_type == 'A' and any(kw in taxonomy for kw in VIRUS_CONFIG['influenza']['types']['A']['taxonomy_keywords'])) or
                    (virus_type == 'B' and any(kw in taxonomy for kw in VIRUS_CONFIG['influenza']['types']['B']['taxonomy_keywords'])) or
                    (virus_type == 'C' and (
                        any(kw in taxonomy for kw in VIRUS_CONFIG['influenza']['types']['A']['taxonomy_keywords']) or
                        any(kw in taxonomy for kw in VIRUS_CONFIG['influenza']['types']['B']['taxonomy_keywords']))
                    )
                ):
                    pass  # Keep processing
                else:
                    continue  # Skip strain if it doesn't match influenza subtype
            else:
                # Match for non-influenza viruses using taxonomy keywords
                if not any(kw in taxonomy for kw in VIRUS_CONFIG[virus_key]['taxonomy_keywords']):
                    continue  # Skip if taxonomy doesn't match

            # --- Handle non-segmented viruses (e.g., HIV, RSV) ---
            if not VIRUS_CONFIG[virus_key].get('segmented', True):
                regions = strain.get('regions', [])
                if regions:
                    region = regions[0]  # Non-segmented viruses usually have one region
                    depth = region.get('depthOfCoverage', 0)
                    coverage = region.get('coveragePercentage', 0)

                    # Check if region passes coverage thresholds
                    passes = depth >= min_depth and coverage >= min_cov

                    results.append({
                        'LimsID': lims_id,
                        'VirusType': virus_key.upper(),
                        'coveragePercentage': strain.get('coveragePercentage'),
                        'depthOfCoverage': strain.get('depthOfCoverage'),
                        'ntIdentity': strain.get('ntIdentity'),
                        'subTypeConclusion': strain.get('subTypeConclusion'),
                        'Submit to GISAID': "Yes" if passes else "No",
                        'Is Control': "Yes" if "PC" in lims_id or "NC" in lims_id or "Neg" in lims_id or "Pos" in lims_id else "No",
                        'ReferenceSequenceId': region.get('referenceSequenceId')
                    })
                continue  # Move to next strain

            # --- Handle segmented viruses (e.g., Influenza A/B) ---
            results.append({
                'LimsID': lims_id,
                'coveragePercentage': strain.get('coveragePercentage'),
                'depthOfCoverage': strain.get('depthOfCoverage'),
                'subTypeConclusion': strain.get('subTypeConclusion'),
                'virusType': 'A' if "alpha" in taxonomy else 'B',
                'regions': [{
                    'segment': r.get('segment'),
                    'depthOfCoverage': r.get('depthOfCoverage'),
                    'coveragePercentage': r.get('coveragePercentage'),
                    'referenceSequenceId': r.get('referenceSequenceId')
                } for r in strain.get('regions', [])]
            })

        return results if results else None

    except Exception as e:
        logging.error(f"Error parsing {json_file}: {e}")
        logging.debug(f"No {virus_type if virus_key == 'influenza' else virus_key} viruses found in {json_file}")
        return None

def evaluate_segment(regions, segment_num, virus_type, min_depth, min_cov):
    """
    Evaluate whether a specific viral genome segment passes defined quality thresholds.

    Parameters:
        regions (list): List of dictionaries, each containing data for a viral genome region.
        segment_num (int): Segment number to evaluate (e.g., 4 for HA or 6 for NA in Influenza).
        virus_type (str): Virus subtype (e.g., 'A' or 'B'). Not directly used here but can aid in debugging.
        min_depth (float): Minimum required depth of coverage.
        min_cov (float): Minimum required percentage coverage.

    Returns:
        str:
            - "Pass" if segment meets both thresholds
            - "Fail (Depth: X, Cov: Y%)" if segment is found but below thresholds
            - "Not found" if the segment isn't present in the regions list
    """
    segment_prefix = f"segment {segment_num}"

    for region in regions:
        segment_name = region.get('segment', '').lower()

        # Match the segment by prefix (e.g., "segment 4")
        if segment_name.startswith(segment_prefix.lower()):
            depth = region.get('depthOfCoverage', 0)
            coverage = region.get('coveragePercentage', 0)

            # Check thresholds
            if depth >= min_depth and coverage >= min_cov:
                return "Pass"
            else:
                return f"Fail (Depth: {depth:.1f}, Cov: {coverage:.1f}%)"

    return "Not found"

def sort_segment_key(segment_name):
    """
    Extract the numeric segment identifier from a segment label for sorting purposes.

    Useful when segment names follow patterns like "Segment 1", "Segment RNA 2", etc.

    Parameters:
        segment_name (str): Name of the segment (e.g., "Segment 4 (HA)").

    Returns:
        int:
            - Extracted segment number as integer if found
            - 999 as fallback to send unrecognized or missing segments to the end of sort
    """
    try:
        # Normalize string and search for a number following "segment" or "segment rna"
        segment_name = segment_name.strip().lower()
        match = re.search(r'segment(?:\s*rna)?\s*(\d+)', segment_name)

        if match:
            return int(match.group(1))  # Return extracted segment number
    except Exception:
        pass

    return 999  # Default fallback value for sorting

def process_influenza_files(temp_folder, influenza_type, min_depth, min_cov):
    """
    Process Genome Detective assignment JSON files for influenza samples 
    (segmented virus) and extract per-segment and per-sample quality metrics.

    For each influenza sample:
        - Segments are sorted numerically
        - Key metrics (depth, coverage, reference) are extracted
        - Segment 4 (HA) and Segment 6 (NA) are used to determine GISAID eligibility

    Parameters:
        temp_folder (str): Path to folder containing `.assignments.json` files.
        influenza_type (str): Subtype selection by user ('A', 'B', or 'C').
        min_depth (float): Minimum coverage depth for segment to pass.
        min_cov (float): Minimum genome coverage percentage for segment to pass.

    Returns:
        tuple:
            - pd.DataFrame: Per-segment information (rows per LIMS ID and segment)
            - pd.DataFrame: GISAID submission summary (one row per LIMS ID)
    """
    segment_table_rows = []  # Stores per-segment metrics
    gisaid_table_rows = []   # Stores per-sample GISAID eligibility summary

    # Iterate over each .assignments.json file
    for json_file in Path(temp_folder).glob('*.assignments.json'):
        influenza_data_list = parse_json_file(json_file, influenza_type, 'influenza', min_depth, min_cov)
        if not influenza_data_list:
            continue  # Skip if no valid influenza strain detected

        for influenza_data in influenza_data_list:
            # Sort segments numerically for consistent output order
            sorted_regions = sorted(
                influenza_data['regions'], 
                key=lambda x: sort_segment_key(x.get('segment', ''))
            )

            # Collect segment metrics
            for region in sorted_regions:
                segment_num = sort_segment_key(region.get('segment', ''))
                gene_name = VIRUS_CONFIG['influenza']['types'][
                    influenza_data['virusType']
                ]['segments'].get(segment_num, f"Segment {segment_num}")
                
                segment_table_rows.append({
                    'LimsID': influenza_data['LimsID'],
                    'VirusType': influenza_data['virusType'],
                    'Segment': f"Segment {segment_num} ({gene_name})",
                    'DepthOfCoverage': region['depthOfCoverage'],
                    'CoveragePercentage': region['coveragePercentage'],
                    'ReferenceSequenceId': region.get('referenceSequenceId')
                })

            # Evaluate Segment 4 (HA) and Segment 6 (NA) quality
            segment4 = evaluate_segment(influenza_data['regions'], 4, influenza_data['virusType'], min_depth, min_cov)
            segment6 = evaluate_segment(influenza_data['regions'], 6, influenza_data['virusType'], min_depth, min_cov)
            submit_to_gisaid = "Yes" if "Pass" in segment4 and "Pass" in segment6 else "No"

            # Check if sample is a control (positive/negative)
            is_control = "Yes" if "PC" in influenza_data['LimsID'] or "NC" in influenza_data['LimsID'] or "Neg" in influenza_data['LimsID'] or "Pos" in influenza_data['LimsID'] else "No"

            # Build GISAID submission summary
            gisaid_table_rows.append({
                'LimsID': influenza_data['LimsID'],
                'VirusType': influenza_data['virusType'],
                f"Segment 4 ({VIRUS_CONFIG['influenza']['types'][influenza_data['virusType']]['segments'][4]})": segment4,
                f"Segment 6 ({VIRUS_CONFIG['influenza']['types'][influenza_data['virusType']]['segments'][6]})": segment6,
                'Submit to GISAID': submit_to_gisaid,
                'Subtype': influenza_data['subTypeConclusion'],
                'Is Control': is_control
            })

    return pd.DataFrame(segment_table_rows), pd.DataFrame(gisaid_table_rows)

def process_non_segmented_virus(temp_folder, virus_key, min_depth, min_cov):
    """
    Process Genome Detective assignment JSON files for non-segmented viruses
    (e.g., HIV, RSV) and collect metrics for GISAID submission.

    Parameters:
        temp_folder (str): Path to folder containing `.assignments.json` files.
        virus_key (str): Virus key used in VIRUS_CONFIG (e.g., 'hiv', 'rsv').
        min_depth (float): Minimum required depth of coverage.
        min_cov (float): Minimum required genome coverage percentage.

    Returns:
        pd.DataFrame: Summary table with one row per sample, including
                      GISAID submission eligibility and key metrics.
    """
    gisaid_table_rows = []

    for json_file in Path(temp_folder).glob('*.assignments.json'):
        virus_data_list = parse_json_file(json_file, '', virus_key, min_depth, min_cov)
        if not virus_data_list:
            continue  # Skip if no matching virus found

        for virus_data in virus_data_list:
            gisaid_table_rows.append({
                'LimsID': virus_data['LimsID'],
                'VirusType': virus_data['VirusType'],
                'CoveragePercentage': virus_data['coveragePercentage'],
                'DepthOfCoverage': virus_data['depthOfCoverage'],
                'Submit to GISAID': virus_data['Submit to GISAID'],
                'Subtype': virus_data['subTypeConclusion'],
                'Is Control': virus_data['Is Control'],
                'ReferenceSequenceId': virus_data['ReferenceSequenceId']
            })

    return pd.DataFrame(gisaid_table_rows)

def extract_genome_info(temp_folder, virus_type, virus_key):
    """
    Extract full-genome summary information for each strain from JSON assignment files.

    For each sample, retrieves:
        - LIMS ID
        - Virus type (A/B/HIV/etc.)
        - Taxonomic name
        - Depth of coverage
        - Genome coverage percentage
        - Number of reads

    Parameters:
        temp_folder (str): Path to folder containing `.assignments.json` files.
        virus_type (str): Subtype for influenza ('A', 'B', or 'C').
        virus_key (str): Lookup key in VIRUS_CONFIG (e.g., 'influenza', 'hiv').

    Returns:
        pd.DataFrame: Full-genome info per matching strain/sample.
    """
    full_genome_data = []

    for json_file in Path(temp_folder).glob('*.assignments.json'):
        try:
            with open(json_file, 'r') as f:
                data = json.load(f)

            lims_id = normalize_limsid(os.path.basename(json_file).split('.')[0])
            strains = data.get('data', {}).get('attributes', {}).get('strains', [])

            for strain in strains:
                taxonomy = strain.get('taxonomyName', '').lower()

                # Match correct virus for influenza
                if virus_key == 'influenza':
                    is_a = any(kw in taxonomy for kw in VIRUS_CONFIG['influenza']['types']['A']['taxonomy_keywords'])
                    is_b = any(kw in taxonomy for kw in VIRUS_CONFIG['influenza']['types']['B']['taxonomy_keywords'])
                    if (virus_type == 'A' and is_a) or \
                       (virus_type == 'B' and is_b) or \
                       (virus_type == 'C' and (is_a or is_b)):
                        full_genome_data.append({
                            'LimsID': lims_id,
                            'VirusType': 'A' if "alpha" in taxonomy else 'B',
                            'taxonomyName': strain.get('taxonomyName'),
                            'numberOfReads': strain.get('numberOfReads'),
                            'depthOfCoverage': strain.get('depthOfCoverage'),
                            'coveragePercentage': strain.get('coveragePercentage'),
                            'ntIdentity': strain.get('ntIdentity')
                        })
                else:
                    # Match non-influenza virus based on keywords
                    if any(kw in taxonomy for kw in VIRUS_CONFIG[virus_key]['taxonomy_keywords']):
                        full_genome_data.append({
                            'LimsID': lims_id,
                            'VirusType': virus_key.upper(),
                            'taxonomyName': strain.get('taxonomyName'),
                            'numberOfReads': strain.get('numberOfReads'),
                            'depthOfCoverage': strain.get('depthOfCoverage'),
                            'coveragePercentage': strain.get('coveragePercentage'),
                            'ntIdentity': strain.get('ntIdentity')
                        })

        except Exception as e:
            logging.error(f"Failed to parse {json_file}: {e}")

    return pd.DataFrame(full_genome_data)

def prompt_for_multiple_labs():
    """
    Prompt the user to indicate whether the dataset contains samples from multiple labs.

    Returns:
        bool: True if multiple labs, False otherwise.
    """
    print("\nAre the samples from multiple labs?")
    choice = input("Enter [y/n]: ").lower().strip()
    while choice not in ['y', 'n']:
        print("Invalid choice. Please enter y or n.")
        choice = input("Enter [y/n]: ").lower().strip()
    return choice == 'y'

def get_lab_info(multiple_labs, all_lims_ids):
    """
    Collect lab assignment information for each sample, either interactively or as a single group.

    Parameters:
        multiple_labs (bool): Whether samples come from multiple labs.
        all_lims_ids (list): List of LIMS IDs that need lab assignment.

    Returns:
        dict: Mapping of lab name to list of LIMS IDs assigned to that lab.
    """
    lab_info = {}
    assigned_lims_ids = set()

    if multiple_labs:
        print("\nYou will now assign samples to different labs.")
        print("At any point, you can:")
        print(" - Enter 'list' to see remaining samples")
        print(" - Enter 'done' to finish assigning")
        print(" - Enter 'labname start-end' to assign a range (e.g., UCT K001-K003)")

        while True:
            remaining = sorted(set(all_lims_ids) - assigned_lims_ids)
            if not remaining:
                print("\nAll samples have been assigned to labs.")
                break

            print(f"\n{len(remaining)} samples remaining")
            entry = input("Enter lab assignment (or 'list'/'done'): ").strip()

            if entry.lower() == 'list':
                print("\nRemaining samples:", ", ".join(remaining))
                continue

            if entry.lower() == 'done':
                print(f"\n{len(remaining)} samples will not be assigned to any lab:")
                print(", ".join(remaining))
                confirm = input("Confirm this is correct? [y/n]: ").lower().strip()
                if confirm == 'y':
                    break
                else:
                    continue

            try:
                # Expecting format like: "UCT K001-K003"
                parts = entry.split()
                if len(parts) < 2:
                    print("Invalid format. Use: LabName StartID-EndID")
                    continue

                lab_name = parts[0]
                range_str = parts[1]
                start, end = range_str.split('-')

                # Use regex to extract letter+number format
                pattern = re.compile(r'(\D+)(\d+)')
                start_match = pattern.match(start)
                end_match = pattern.match(end)

                if not start_match or not end_match:
                    print("Invalid ID format. Use format like K001-K003")
                    continue

                prefix = start_match.group(1)
                start_num = int(start_match.group(2))
                end_num = int(end_match.group(2))

                # Match samples in that numeric range with the same prefix
                matched_ids = []
                for lims_id in remaining:
                    lims_match = pattern.match(lims_id)
                    if lims_match and lims_match.group(1) == prefix:
                        num = int(lims_match.group(2))
                        if start_num <= num <= end_num:
                            matched_ids.append(lims_id)

                if not matched_ids:
                    print(f"No samples found in range {start}-{end}")
                    continue

                print(f"\nAssigning {len(matched_ids)} samples to {lab_name}:")
                print(", ".join(matched_ids))
                confirm = input("Confirm this assignment? [y/n]: ").lower().strip()
                if confirm != 'y':
                    continue

                lab_info[lab_name] = matched_ids
                assigned_lims_ids.update(matched_ids)

            except Exception as e:
                print(f"Error processing input: {e}")
                continue

    else:
        # For single-lab assignment, offer common lab options
        print("\nCommon lab names:")
        print("1. NHLS-SU")
        print("2. NHLS-UCT")
        print("3. Pathcare (PATH)")
        print("4. Other")
        choice = input("Select lab (1-4): ").strip()

        if choice == '1':
            lab_name = 'NHLS-SU'
        elif choice == '2':
            lab_name = 'NHLS-UCT'
        elif choice == '3':
            lab_name = 'PATH'
        else:
            lab_name = input("Enter lab name: ").strip()

        lab_info[lab_name] = all_lims_ids

    return lab_info

def prompt_user_for_header_format(virus_type, virus_key, all_lims_ids, lab_info=None):
    """
    Prompt the user to confirm or customize the FASTA header format.

    Header strings use placeholders like:
        - <LimsID>: Sample identifier
        - <lab>: Assigned lab name
        - <gene>: Segment gene (for influenza)
        - <isolate>: Virus subtype (A or B)

    Parameters:
        virus_type (str): Type of virus ('A', 'B', or 'C' for influenza).
        virus_key (str): Key in the VIRUS_CONFIG dictionary.
        all_lims_ids (list): List of LIMS IDs eligible for submission.
        lab_info (dict, optional): Mapping of lab names to LIMS IDs.

    Returns:
        tuple:
            - str: Header template string for FASTA output.
            - dict: Lab info mapping (lab name → list of LIMS IDs).
    """
    # If no lab info provided, prompt the user interactively
    if lab_info is None:
        multiple_labs = prompt_for_multiple_labs()
        lab_info = get_lab_info(multiple_labs, all_lims_ids)

    # Define default FASTA header template
    if virus_key == 'influenza':
        if virus_type == 'C':
            # For combined A+B, isolate subtype dynamically during writing
            default_header = "><isolate>/South Africa/<lab>-CERI-<LimsID>/2025_<gene>"
        else:
            default_header = f">{virus_type}/South Africa/<lab>-CERI-<LimsID>/2025_<gene>"
    else:
        default_header = f">{VIRUS_CONFIG[virus_key]['common_name']}/South Africa/<lab>-CERI-<LimsID>/2025"

    # Display and confirm header format
    print(f"\nDefault header format: {default_header}")
    confirm = input("Are you happy with this format? (y/n): ").lower()

    if confirm != 'y':
        # Let user define a custom format depending on virus type
        if virus_key == 'influenza':
            custom_format = input("Enter custom format using <LimsID>, <gene>, <lab> and <isolate>: ")
        else:
            custom_format = input("Enter custom format using <LimsID> and <lab>: ")
        return custom_format, lab_info
    else:
        return default_header, lab_info

def get_lab_for_limsid(limsid, lab_info):
    """
    Retrieve the lab name corresponding to a specific LIMS ID.

    Parameters:
        limsid (str): The sample LIMS ID.
        lab_info (dict): Dictionary mapping lab names to lists of LIMS IDs.

    Returns:
        str: The lab name corresponding to the LIMS ID, or "UNKNOWN" if not found.
    """
    if len(lab_info) == 1:
        # If there's only one lab, return it directly
        return next(iter(lab_info.keys()))

    # Otherwise, search for which lab contains the LIMS ID
    for lab_name, lims_ids in lab_info.items():
        if limsid in lims_ids:
            return lab_name

    return "UNKNOWN"

def get_segment_number_and_gene(header_line, virus_name, virus_type):
    """
    Extract segment number and corresponding gene name from a segment label.

    Expected formats:
        - "Segment 4"
        - "Segment 4 (HA)"
        - "4" (as fallback)

    Parameters:
        header_line (str): Segment label or description.
        virus_name (str): Key in VIRUS_CONFIG (e.g., 'influenza').
        virus_type (str): Subtype (e.g., 'A' or 'B').

    Returns:
        tuple:
            - int or None: Segment number
            - str or None: Gene name (e.g., 'HA', 'NA', etc.)
    """
    try:
        # Match labels like "Segment 4" or "Segment RNA 4"
        match = re.search(r'[Ss]egment\s*(\d+)', str(header_line))
        if match:
            segment_num = int(match.group(1))
            gene_name = VIRUS_CONFIG[virus_name]['types'][virus_type]['segments'].get(
                segment_num, f"Segment{segment_num}"
            )
            return segment_num, gene_name

        # Fallback: header is just a digit (e.g., "4")
        elif str(header_line).strip().isdigit():
            segment_num = int(header_line.strip())
            gene_name = VIRUS_CONFIG[virus_name]['types'][virus_type]['segments'].get(
                segment_num, f"Segment{segment_num}"
            )
            return segment_num, gene_name

    except Exception as e:
        logging.warning(f"Error parsing segment from header '{header_line}': {e}")

    return None, None

def write_segment_fastas(input_folder, segment_df, gisaid_df, virus_type):
    """
    Write individual FASTA files for each segment of influenza virus samples that passed GISAID thresholds.

    For each eligible sample:
        - Retrieves reference FASTA from alignment output in sample-specific folder
        - Extracts the consensus sequence (always the second sequence in file)
        - Writes one segment per FASTA file with a standardized or user-defined header

    Parameters:
        input_folder (str): Path to the folder containing raw alignment FASTA files.
        segment_df (pd.DataFrame): DataFrame containing segment-level metadata (depth, coverage, etc.).
        gisaid_df (pd.DataFrame): DataFrame containing GISAID submission eligibility and subtype.
        virus_type (str): User-selected virus type ('A', 'B', or 'C' for both A and B).
    """
    eligible = gisaid_df[gisaid_df['Submit to GISAID'] == 'Yes']
    all_lims_ids = eligible['LimsID'].tolist()

    temp2_folder = os.path.join(os.getcwd(), 'Outputs')
    os.makedirs(temp2_folder, exist_ok=True)

    header_template, lab_info = prompt_user_for_header_format(virus_type, 'influenza', all_lims_ids)

    for _, row in eligible.iterrows():
        limsid = row['LimsID']
        original_limsid = limsid

        if "PC" in limsid or "NC" in limsid or "Neg" in limsid or "Pos" in limsid:
            pass
        elif "-" in limsid:
            limsid = limsid.split("-")[0]

        current_virus_type = row['VirusType']
        segments = segment_df[(segment_df['LimsID'] == limsid) & (segment_df['VirusType'] == current_virus_type)]

        if virus_type == 'C':
            dst_folder = os.path.join(temp2_folder, limsid, f'Influenza_{current_virus_type}')
        else:
            dst_folder = os.path.join(temp2_folder, limsid)
        os.makedirs(dst_folder, exist_ok=True)

        sample_input_folder = None
        for folder in os.listdir(input_folder):
            if original_limsid in folder:
                sample_input_folder = os.path.join(input_folder, folder)
                break

        if not sample_input_folder or not os.path.exists(sample_input_folder):
            logging.warning(f"No folder containing {original_limsid} found in {input_folder}")
            continue

        for _, seg_row in segments.iterrows():
            ref_id = seg_row['ReferenceSequenceId']
            depth_id = seg_row['DepthOfCoverage']
            coverage_id = seg_row['CoveragePercentage']
            seg_name = seg_row['Segment']
            seg_num, gene_name = get_segment_number_and_gene(seg_name, 'influenza', current_virus_type)

            if not seg_num:
                logging.warning(f"Could not extract segment number from {seg_name} for sample {limsid}")
                continue

            if depth_id > 10 and coverage_id > 80:
                found = False  # Move declaration outside the file loop
                for file in os.listdir(sample_input_folder):
                    if ref_id in file and file.endswith('alignment-nt.fasta'):
                        src = os.path.join(sample_input_folder, file)
                        dst_file = os.path.join(dst_folder, f'segment{seg_num}.fasta')
                        found = True

                        try:
                            with open(src, 'r') as f:
                                lines = [line.strip() for line in f if line.strip()]

                            headers = [i for i, line in enumerate(lines) if line.startswith('>')]
                            if len(headers) != 2:
                                logging.warning(f"Expected exactly 2 sequences in {src}, found {len(headers)}. Skipping.")
                                continue

                            consensus_start = headers[1]
                            full_sequence = ''.join(lines[consensus_start + 1:])

                            lab_name = get_lab_for_limsid(limsid, lab_info)

                            new_header = header_template \
                                .replace('<LimsID>', limsid) \
                                .replace('<gene>', gene_name) \
                                .replace('<lab>', lab_name) \
                                .replace('<isolate>', f'{current_virus_type}')

                            if not new_header.startswith('>'):
                                new_header = '>' + new_header

                            new_header = new_header.replace(f'>{current_virus_type}/{current_virus_type}/', f'>{current_virus_type}/')
                            new_header = new_header.replace('>AB/', f'>{current_virus_type}/')

                            with open(dst_file, 'w') as out:
                                out.write(f"{new_header}\n{full_sequence}\n")

                            logging.info(f"Written segment fasta: {dst_file}")
                            break

                        except Exception as e:
                            logging.error(f"Failed to process {src}: {e}")
                            continue

                if not found:
                    logging.warning(f"No alignment file found for reference ID {ref_id} in {sample_input_folder}")
       

def write_virus_fastas(input_folder, df, virus_key, lab_info):
    """
    Write FASTA files for non-segmented virus samples that passed GISAID criteria.

    For each sample:
        - Locate the corresponding alignment file in the sample-specific folder
        - Extract the consensus sequence
        - Write it with a user-defined or default FASTA header

    Parameters:
        input_folder (str): Directory where the alignment FASTA files are stored.
        df (pd.DataFrame): GISAID summary dataframe for the virus.
        virus_key (str): Key from VIRUS_CONFIG (e.g., 'hiv', 'rsv').
        lab_info (dict): Mapping of lab names to LIMS IDs.
    """
    temp_folder = os.path.join(os.getcwd(), 'Outputs')
    os.makedirs(temp_folder, exist_ok=True)

    # Display default header format for confirmation
    print(f"\nDefault header format for {VIRUS_CONFIG[virus_key]['name']}:")
    default_header = f">{VIRUS_CONFIG[virus_key]['common_name']}/South Africa/<lab>-CERI-<LimsID>/2025"
    print(default_header)

    # Allow custom header override
    confirm = input("Use this format? (y/n): ").lower()
    if confirm != 'y':
        header_template = input("Enter custom format using <LimsID> and <lab>: ")
    else:
        header_template = default_header

    # Process only GISAID-eligible samples
    eligible_samples = df[df['Submit to GISAID'] == 'Yes']

    for _, row in eligible_samples.iterrows():
        limsid = row['LimsID']
        ref_id = row['ReferenceSequenceId']

        # Create sample output folder
        sample_folder = os.path.join(temp_folder, limsid)
        os.makedirs(sample_folder, exist_ok=True)

        sample_input_folder = None
        for folder in os.listdir(input_folder):
            if limsid in folder:
                sample_input_folder = os.path.join(input_folder, folder)
                break

        if not sample_input_folder or not os.path.exists(sample_input_folder):
            logging.warning(f"No folder containing {limsid} found in {input_folder}")
            continue

        found = False
        # Search only in the found sample folder
        for file in os.listdir(sample_input_folder):
            if ref_id in file and file.endswith('alignment-nt.fasta'):
                src = os.path.join(sample_input_folder, file)
                dst_file = os.path.join(sample_folder, f"{limsid}_{virus_key}.fasta")
                found = True

                try:
                    # Extract sequence after second FASTA header
                    with open(src, 'r') as f:
                        lines = [line.strip() for line in f if line.strip()]

                    headers = [i for i, line in enumerate(lines) if line.startswith('>')]
                    if len(headers) < 2:
                        logging.warning(f"Expected 2 headers in {src}, found {len(headers)}. Skipping.")
                        continue

                    sequence = ''.join(lines[headers[1] + 1:])

                    # Replace template values in header
                    lab_name = get_lab_for_limsid(limsid, lab_info)
                    new_header = header_template\
                        .replace('<LimsID>', limsid)\
                        .replace('<lab>', lab_name)

                    # Write final FASTA file
                    with open(dst_file, 'w') as out:
                        out.write(f"{new_header}\n{sequence}\n")

                    logging.info(f"Written FASTA file: {dst_file}")
                    break  # Found the file, no need to keep searching

                except Exception as e:
                    logging.error(f"Failed to process {src}: {e}")
                    continue

        if not found:
            logging.warning(f"No alignment file found for reference ID {ref_id} in {sample_input_folder}")

def concatenate_segments_per_sample(virus_type):
    """
    Concatenate individual influenza segment FASTA files per sample into one multi-segment file.

    Handles both:
        - Single-subtype samples (A or B)
        - Dual-subtype setup ('C'/‘AB’) by merging `Influenza_A` and `Influenza_B` subfolders separately.

    Output is written to: 
        Outputs/<LimsID>[/Influenza_A or B]/<LimsID>_all_segments.fasta

    Parameters:
        virus_type (str): 'A', 'B', or 'AB' (for both A and B as used in "C" mode)
    """
    temp2_folder = os.path.join(os.getcwd(), 'Outputs')

    for full_folder_name in os.listdir(temp2_folder):
        # Normalize folder name by stripping suffix if needed
        if "PC" in full_folder_name or "NC" in full_folder_name or "Neg" in full_folder_name or "Pos" in full_folder_name:
            limsid_folder = full_folder_name
        else:
            limsid_folder = full_folder_name.split('-')[0]

        sample_folder = os.path.join(temp2_folder, limsid_folder)
        if not os.path.isdir(sample_folder):
            continue

        if virus_type == 'AB':
            # Process both Influenza_A and Influenza_B subfolders
            for influ_type in ['A', 'B']:
                type_folder = os.path.join(sample_folder, f'Influenza_{influ_type}')
                if not os.path.exists(type_folder):
                    continue

                segment_files = [
                    os.path.join(type_folder, f)
                    for f in os.listdir(type_folder)
                    if f.startswith('segment') and f.endswith('.fasta')
                ]
                segment_files.sort()

                if segment_files:
                    concatenated_fasta = os.path.join(type_folder, f"{limsid_folder}_all_segments.fasta")
                    with open(concatenated_fasta, 'w') as outfile:
                        for seg_file in segment_files:
                            with open(seg_file, 'r') as infile:
                                content = infile.read().strip()
                                outfile.write(content + '\n')
                    logging.info(f"Concatenated segments for sample {limsid_folder} (type {influ_type}) into {concatenated_fasta}")

        else:
            # Single-type influenza sample
            segment_files = [
                os.path.join(sample_folder, f)
                for f in os.listdir(sample_folder)
                if f.startswith('segment') and f.endswith('.fasta')
            ]
            segment_files.sort()

            if segment_files:
                concatenated_fasta = os.path.join(sample_folder, f"{limsid_folder}_all_segments.fasta")
                with open(concatenated_fasta, 'w') as outfile:
                    for seg_file in segment_files:
                        with open(seg_file, 'r') as infile:
                            content = infile.read().strip()
                            outfile.write(content + '\n')
                logging.info(f"Concatenated segments for sample {limsid_folder} into {concatenated_fasta}")

def concatenate_all_samples_fasta(virus_type, output_filename_prefix="all_samples_combined"):
    """
    Concatenate all per-sample full-genome FASTA files into a single master file.

    For influenza:
        - If virus_type is 'AB', generates two files:
            - all_samples_combined_A.fasta
            - all_samples_combined_B.fasta

    For other viruses:
        - One combined file is created: all_samples_combined.fasta

    Parameters:
        virus_type (str): 'A', 'B', or 'AB' for influenza, or other virus type keys.
        output_filename_prefix (str): Prefix for output FASTA file(s). Default is 'all_samples_combined'.
    """
    temp2_folder = os.path.join(os.getcwd(), 'Outputs')

    if virus_type == 'AB':
        # Process A and B types separately
        for influ_type in ['A', 'B']:
            combined_fasta_path = os.path.join(temp2_folder, f"{output_filename_prefix}_{influ_type}.fasta")
            with open(combined_fasta_path, 'w') as outfile:
                for limsid_folder in os.listdir(temp2_folder):
                    sample_folder = os.path.join(temp2_folder, limsid_folder)
                    if not os.path.isdir(sample_folder):
                        continue

                    if "PC" in limsid_folder or "NC" in limsid_folder or "Neg" in limsid_folder or "Pos" in limsid_folder:
                        logging.info(f"Skipping control sample: {limsid_folder}")
                        continue

                    type_folder = os.path.join(sample_folder, f'Influenza_{influ_type}')
                    if not os.path.exists(type_folder):
                        continue

                    all_segments_fasta = os.path.join(type_folder, f"{limsid_folder}_all_segments.fasta")
                    if not os.path.exists(all_segments_fasta):
                        logging.warning(f"Missing concatenated fasta for sample {limsid_folder}: {all_segments_fasta}")
                        continue

                    with open(all_segments_fasta, 'r') as infile:
                        content = infile.read().strip()
                        outfile.write(f">Sample_{limsid_folder}\n")
                        outfile.write(content + '\n')

            logging.info(f"Combined fasta for all {influ_type} samples written to {combined_fasta_path}")

    else:
        # Combine all into a single FASTA
        combined_fasta_path = os.path.join(temp2_folder, f"{output_filename_prefix}.fasta")
        with open(combined_fasta_path, 'w') as outfile:
            for limsid_folder in os.listdir(temp2_folder):
                sample_folder = os.path.join(temp2_folder, limsid_folder)
                if not os.path.isdir(sample_folder):
                    continue

                if "PC" in limsid_folder or "NC" in limsid_folder or "Neg" in limsid_folder or "Pos" in limsid_folder:
                    logging.info(f"Skipping control sample: {limsid_folder}")
                    continue

                fasta_pattern = os.path.join(sample_folder, f"{limsid_folder}*.fasta")
                matching_fastas = glob.glob(fasta_pattern)

                if not matching_fastas:
                    logging.warning(f"Missing concatenated fasta for sample {limsid_folder}: No match for pattern {fasta_pattern}")
                    continue

                # Use the first match (or loop over all if multiple are expected)
                fasta_file_path = matching_fastas[0]
                with open(fasta_file_path, 'r') as infile:
                    content = infile.read().strip()
                    outfile.write(content + '\n')


        logging.info(f"Combined fasta for all samples written to {combined_fasta_path}")

def main(input_folder):
    """
    Main entry point for the viral analysis pipeline.

    This function:
    1. Prompts the user to select virus type and coverage thresholds.
    2. Extracts and filters data from Genome Detective assignment files.
    3. Generates summary statistics and CSV reports.
    4. Writes per-sample FASTA files and concatenates segments for influenza.
    5. Handles both segmented (influenza A/B) and non-segmented viruses.

    Parameters:
        input_folder (str): Path to the folder containing JSON assignment files and alignment FASTAs.

    Returns:
        bool: True on successful completion.
    """
    logging.info("Starting analysis")

    # --- Step 1: Collect user input ---
    virus_type, virus_key = get_virus_type()
    virus_type_header = (
        'A' if virus_type == 'A'
        else 'B' if virus_type == 'B'
        else 'AB' if virus_key == 'influenza'
        else virus_key
    )
    min_depth, min_cov = get_coverage_thresholds()

    # --- Step 2: Preprocess input files ---
    print("STEP 2")
    temp_folder = create_temp_folder(input_folder)

    # Parse JSON files depending on virus type
    if virus_key == 'influenza':
        segment_df, gisaid_df = process_influenza_files(temp_folder, virus_type, min_depth, min_cov)
    else:
        gisaid_df = process_non_segmented_virus(temp_folder, virus_key, min_depth, min_cov)
        segment_df = None

    # Full genome info summary (used for all viruses)
    full_genome_df = extract_genome_info(temp_folder, virus_type, virus_key)

    # --- Step 3: Validation check for detected viruses ---
    if virus_key == 'influenza':
        if virus_type in ['A', 'B'] and len(gisaid_df) == 0:
            logging.error(f"Your chosen virus type (Influenza {virus_type}) was not detected in any samples")
            exit(1)
        elif virus_type == 'C' and (
            len(gisaid_df[gisaid_df['VirusType'] == 'A']) == 0 and
            len(gisaid_df[gisaid_df['VirusType'] == 'B']) == 0
        ):
            logging.error("Neither Influenza A nor B were detected in any samples")
            exit(1)
    else:
        if len(gisaid_df) == 0:
            logging.error(f"Your chosen virus ({VIRUS_CONFIG[virus_key]['name']}) was not detected in any samples")
            exit(1)

    # --- Step 4: Generate summary CSVs ---
    output_dir = os.path.join(os.getcwd(), 'Summary_files')
    os.makedirs(output_dir, exist_ok=True)

    if virus_key == 'influenza':
        if virus_type == 'C':  # User selected both A and B
            for v_type in ['A', 'B']:
                # Filter dataframes per subtype
                type_segment_df = segment_df[segment_df['VirusType'] == v_type].sort_values(by='LimsID')
                type_gisaid_df = gisaid_df[gisaid_df['VirusType'] == v_type].sort_values(by='LimsID')
                type_genome_df = full_genome_df[full_genome_df['VirusType'] == v_type].sort_values(by='LimsID')

                # Write output CSVs
                type_segment_df.to_csv(os.path.join(output_dir, f'influenza_{v_type.lower()}_segments.csv'), index=False)
                type_gisaid_df.to_csv(os.path.join(output_dir, f'gisaid_submission_status_{v_type.lower()}.csv'), index=False)
                type_genome_df.to_csv(os.path.join(output_dir, f'full_genome_info_{v_type.lower()}.csv'), index=False)

                # Log stats
                eligible_count = len(type_gisaid_df[type_gisaid_df['Submit to GISAID'] == 'Yes'])
                control_count = len(type_gisaid_df[(type_gisaid_df['Submit to GISAID'] == 'Yes') &
                                                   (type_gisaid_df['Is Control'] == 'Yes')])
                logging.info(f"Influenza {v_type} - Total samples: {len(type_gisaid_df)}, "
                             f"GISAID-eligible: {eligible_count} (including {control_count} controls)")
        else:
            # A or B only
            segment_df = segment_df.sort_values(by='LimsID')
            gisaid_df = gisaid_df.sort_values(by='LimsID')
            full_genome_df = full_genome_df.sort_values(by='LimsID')

            segment_df.to_csv(os.path.join(output_dir, f'influenza_{virus_type.lower()}_segments.csv'), index=False)
            gisaid_df.to_csv(os.path.join(output_dir, f'gisaid_submission_status_{virus_type.lower()}.csv'), index=False)
            full_genome_df.to_csv(os.path.join(output_dir, f'full_genome_info_{virus_type.lower()}.csv'), index=False)

            eligible_count = len(gisaid_df[gisaid_df['Submit to GISAID'] == 'Yes'])
            control_count = len(gisaid_df[(gisaid_df['Submit to GISAID'] == 'Yes') &
                                          (gisaid_df['Is Control'] == 'Yes')])
            logging.info(f"Total samples: {len(gisaid_df)}, GISAID-eligible: {eligible_count} (including {control_count} controls)")
    else:
        # Non-influenza viruses
        gisaid_df = gisaid_df.sort_values(by='LimsID')
        full_genome_df = full_genome_df.sort_values(by='LimsID')

        gisaid_df.to_csv(os.path.join(output_dir, f'gisaid_submission_status_{virus_key}.csv'), index=False)
        full_genome_df.to_csv(os.path.join(output_dir, f'full_genome_info_{virus_key}.csv'), index=False)

        eligible_count = len(gisaid_df[gisaid_df['Submit to GISAID'] == 'Yes'])
        control_count = len(gisaid_df[(gisaid_df['Submit to GISAID'] == 'Yes') &
                                      (gisaid_df['Is Control'] == 'Yes')])
        logging.info(f"Total samples: {len(gisaid_df)}, GISAID-eligible: {eligible_count} (including {control_count} controls)")

    # --- Step 5: Write FASTA files and concatenate outputs ---
    if virus_key == 'influenza':
        write_segment_fastas(input_folder, segment_df, gisaid_df, virus_type)
        concatenate_segments_per_sample(virus_type_header)
        concatenate_all_samples_fasta(virus_type_header)
    else:
        # Prompt for lab info and write FASTAs
        all_lims_ids = gisaid_df['LimsID'].tolist()
        multiple_labs = prompt_for_multiple_labs()
        lab_info = get_lab_info(multiple_labs, all_lims_ids)
        write_virus_fastas(input_folder, gisaid_df, virus_key, lab_info)
        concatenate_all_samples_fasta(virus_type_header)

    return True

if __name__ == "__main__":
    import argparse
    import logging

    # Set up command-line argument parsing
    parser = argparse.ArgumentParser(
        description=' Viral GISAID Submission Pipeline\n'
                    'Processes Genome Detective outputs and prepares submission files.'
    )
    parser.add_argument(
        'input_folder',
        help='Path to input folder containing LimsID-named subdirectories with .assignments.json and alignment FASTAs'
    )
    args = parser.parse_args()

    # Optional: Set up logging (if not already configured globally)
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s [%(levelname)s] %(message)s',
        handlers=[logging.StreamHandler()]
    )

    # Run the main pipeline
    success = main(args.input_folder)

    if not success:
        logging.error("Pipeline failed.")
        exit(1)
    else:
        logging.info("Pipeline completed successfully.")

