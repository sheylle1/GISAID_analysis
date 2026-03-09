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
- Handles both old and new Genome Detective JSON formats
- Prevents duplicate processing of samples
- Extracts true unique ID (C0XXXX or K0XXXX) from sample identifiers
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
        'taxonomy_keywords': ['orthopneumovirus hominis', 'human respiratory syncytial virus'],
        'common_name': 'RSV'
    },
    'covid': {
        'name': 'SARS-CoV-2',
        'segmented': False,
        'taxonomy_keywords': [
            'severe acute respiratory syndrome-related coronavirus',  # Old format
            'betacoronavirus pandemicum'  # New format
        ],
        'common_name': 'hCoV-19'
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

def extract_unique_id(sample_id):
    """
    Extract the true unique ID from a sample identifier.
    The unique ID is always a C0XXXX or K0XXXX pattern (C or K followed by numbers).
    
    Parameters:
        sample_id (str): The raw sample ID that may contain prefixes/suffixes
        
    Returns:
        str: The extracted unique ID, or the original if no pattern found
    """
    if not sample_id:
        return sample_id
        
    # Pattern for C0XXXX or K0XXXX (C or K followed by numbers)
    # This matches patterns like C01234, K01234, C001234, etc.
    pattern = r'([CK]0\d+)'
    match = re.search(pattern, sample_id)
    
    if match:
        return match.group(1)
    
    # If no pattern found, return original
    return sample_id

def normalize_limsid(raw_id):
    """
    Normalize a raw LIMS ID by removing suffixes for non-control samples.
    
    Parameters:
        raw_id (str): The raw sample ID from the filename.

    Returns:
        str: Normalized LIMS ID.
    """
    # First, check if it's a control
    if "PC" in raw_id or "NC" in raw_id or "Neg" in raw_id or "Pos" in raw_id or "ERCC" in raw_id:
        return raw_id  # Keep controls as-is
    
    # Check for duplicates with -D suffix
    if "-D" in raw_id:
        return raw_id  # Keep duplicates as-is as they're intentionally separate
    
    # Extract the unique C0 or K0 ID
    unique_id = extract_unique_id(raw_id)
    
    return unique_id

def create_temp_folder(input_folder):
    """Create temporary folder and copy assignment files."""
    temp_folder = os.path.join(os.getcwd(), 'GD_Assignments')
    os.makedirs(temp_folder, exist_ok=True)
    logging.info(f"Created temp folder: {temp_folder}")
    
    # Track unique files to avoid duplicates
    copied_files = set()
    
    for root, _, files in os.walk(input_folder):
        for file in files:
            if file.endswith('.assignments.json'):
                src_path = os.path.join(root, file)
                dest_path = os.path.join(temp_folder, file)
                
                # Skip if we've already copied this file
                if file in copied_files:
                    logging.debug(f"Skipping duplicate file: {file}")
                    continue
                    
                if os.path.abspath(src_path) != os.path.abspath(dest_path):
                    try:
                        with open(src_path, 'rb') as src_file:
                            data = src_file.read()
                            if not data:
                                continue
                            else:
                                with open(dest_path, 'wb') as dest_file:
                                    dest_file.write(data)
                                copied_files.add(file)
                                #logging.info(f"Copied {file} to temp folder")
                    except Exception as e:
                        logging.error(f"Failed to copy {src_path}: {e}")
    return temp_folder

def get_virus_type():
    """
    Prompt the user to select the virus type from a predefined list.
    
    Returns:
        tuple: (virus_type (str), virus_key (str))
    """
    # Display menu of available viruses
    print("\nSelect virus type:")
    virus_list = list(VIRUS_CONFIG.keys())
    for i, (virus_key, virus_data) in enumerate(VIRUS_CONFIG.items(), 1):
        print(f"{i}. {virus_data['name']}")
    
    # Get and validate user input
    choice = input(f"Enter your choice (1-{len(VIRUS_CONFIG)}): ").strip()
    while not choice.isdigit() or int(choice) not in range(1, len(VIRUS_CONFIG) + 1):
        print("Invalid choice. Please try again.")
        choice = input(f"Enter your choice (1-{len(VIRUS_CONFIG)}): ").strip()
    
    # Resolve virus key from numeric selection
    virus_key = virus_list[int(choice) - 1]
    
    # If Influenza, follow up with type-specific prompt (A/B/C)
    if virus_key == 'influenza':
        return get_influenza_type(), virus_key
    
    return virus_key.upper(), virus_key

def get_influenza_type():
    """
    Prompt the user to select which type of influenza virus samples are being analyzed.

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

    Returns:
        tuple: (float, float) minimum depth and coverage percentage
    """
    print("\nSet coverage thresholds for GISAID submission:")
    min_depth = input("Minimum depth of coverage (default 10): ").strip()
    min_cov = input("Minimum coverage percentage (default 80): ").strip()
    
    try:
        min_depth = float(min_depth) if min_depth else 10.0
        min_cov = float(min_cov) if min_cov else 80.0
    except ValueError:
        logging.warning("Invalid input. Using default thresholds (depth=10, cov=80%)")
        min_depth, min_cov = 10.0, 80.0
    
    return min_depth, min_cov

def extract_segment_number_and_gene(segment_str):
    """
    Extract segment number and gene name from segment string.
    Handles both old format ("segment X-NAME") and new format ("X-NAME").
    
    Parameters:
        segment_str (str): Segment string (e.g., "segment 4-HA" or "4-HA")
    
    Returns:
        tuple: (segment_number, gene_name) or (None, None) if parsing fails
    """
    if not segment_str or not isinstance(segment_str, str):
        return None, None
    
    # Remove "segment " prefix if present (for old format)
    cleaned = re.sub(r'^segment\s+', '', segment_str.lower(), flags=re.IGNORECASE)
    
    # Try to match pattern like "4-ha" or "4" (number followed by optional hyphen and gene)
    match = re.match(r'^(\d+)(?:-([a-z0-9]+))?$', cleaned)
    if match:
        segment_num = int(match.group(1))
        gene_name = match.group(2).upper() if match.group(2) else None
        return segment_num, gene_name
    
    return None, None

def parse_json_file(json_file, virus_type, virus_key, min_depth, min_cov):
    """
    Parse a Genome Detective JSON `.assignments.json` file and extract viral information.
    Handles both old and new JSON formats.
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

        # Extract the unique ID from the filename
        base_name = os.path.basename(json_file).replace('.assignments.json', '')
        lims_id = normalize_limsid(base_name)

        # Get list of all strains from the file
        strains = data.get('data', {}).get('attributes', {}).get('strains', [])

        results = []

        for strain in strains:
            taxonomy = strain.get('taxonomyName', '').lower()

            # --- Match strain to correct virus type ---
            if virus_key == 'influenza':
                # Influenza-specific matching based on subtype
                is_a = any(kw in taxonomy for kw in VIRUS_CONFIG['influenza']['types']['A']['taxonomy_keywords'])
                is_b = any(kw in taxonomy for kw in VIRUS_CONFIG['influenza']['types']['B']['taxonomy_keywords'])
                
                if (virus_type == 'A' and is_a) or \
                   (virus_type == 'B' and is_b) or \
                   (virus_type == 'C' and (is_a or is_b)):
                    pass  # Keep processing
                else:
                    continue  # Skip strain if it doesn't match influenza subtype
            else:
                # Match for non-influenza viruses using taxonomy keywords
                matched = False
                for kw in VIRUS_CONFIG[virus_key]['taxonomy_keywords']:
                    if kw.lower() in taxonomy:
                        matched = True
                        break
                
                if not matched:
                    continue  # Skip if taxonomy doesn't match

            # --- Handle non-segmented viruses (e.g., HIV, RSV, COVID) ---
            if not VIRUS_CONFIG[virus_key].get('segmented', True):
                regions = strain.get('regions', [])
                if regions:
                    # For non-segmented viruses, take the best region (highest coverage)
                    best_region = None
                    best_coverage = 0
                    
                    for region in regions:
                        coverage = region.get('coveragePercentage', 0)
                        if coverage > best_coverage:
                            best_coverage = coverage
                            best_region = region
                    
                    if best_region:
                        depth = best_region.get('depthOfCoverage', 0)
                        coverage = best_region.get('coveragePercentage', 0)

                        # Check if region passes coverage thresholds
                        passes = depth >= min_depth and coverage >= min_cov

                        results.append({
                            'LimsID': lims_id,
                            'VirusType': virus_key.upper(),
                            'coveragePercentage': coverage,
                            'depthOfCoverage': depth,
                            'ntIdentity': strain.get('ntIdentity'),
                            'subTypeConclusion': strain.get('subTypeConclusion'),
                            'Submit to GISAID': "Yes" if passes else "No",
                            'Is Control': "Yes" if any(ctrl in lims_id for ctrl in ["PC", "NC", "Neg", "Pos", "ERCC"]) else "No",
                            'ReferenceSequenceId': best_region.get('referenceSequenceId')
                        })
                continue  # Move to next strain

            # --- Handle segmented viruses (e.g., Influenza A/B) ---
            # Determine virus type from taxonomy
            # Determine virus type using config keywords (robust)
            current_virus_type = None

            for vtype, vconfig in VIRUS_CONFIG['influenza']['types'].items():
                for kw in vconfig['taxonomy_keywords']:
                    if kw in taxonomy:
                        current_virus_type = vtype
                        break
                if current_virus_type:
                    break

            if not current_virus_type:
                logging.warning(f"Could not determine virus type from taxonomy: {taxonomy}")
                continue
            
            # Process all regions
            region_data = []
            
            for region in strain.get('regions', []):
                segment_str = region.get('segment', '')
                segment_num, gene = extract_segment_number_and_gene(segment_str)
                
                region_data.append({
                    'segment': segment_str,
                    'segment_num': segment_num,
                    'gene': gene,
                    'depthOfCoverage': region.get('depthOfCoverage', 0),
                    'coveragePercentage': region.get('coveragePercentage', 0),
                    'referenceSequenceId': region.get('referenceSequenceId')
                })

            # Only add if we have regions
            if region_data:
                results.append({
                    'LimsID': lims_id,
                    'coveragePercentage': strain.get('coveragePercentage'),
                    'depthOfCoverage': strain.get('depthOfCoverage'),
                    'subTypeConclusion': strain.get('subTypeConclusion'),
                    'virusType': current_virus_type,
                    'regions': region_data
                })

    

        return results if results else None

    except Exception as e:
        logging.error(f"Error parsing {json_file}: {e}")
        return None

def evaluate_segment(regions, segment_num, virus_type, min_depth, min_cov):
    """
    Evaluate whether a specific viral genome segment passes defined quality thresholds.
    Handles both old and new segment naming formats.
    """
    for region in regions:
        segment_str = region.get('segment', '')
        segment_num_from_str, _ = extract_segment_number_and_gene(segment_str)
        
        if segment_num_from_str == segment_num:
            depth = region.get('depthOfCoverage', 0)
            coverage = region.get('coveragePercentage', 0)

            if depth >= min_depth and coverage >= min_cov:
                return "Pass"
            else:
                return f"Fail (Depth: {depth:.1f}, Cov: {coverage:.1f}%)"

    return "Not found"

def sort_segment_key(region):
    """
    Extract the numeric segment identifier from a region for sorting purposes.
    Handles both old and new formats.
    """
    segment_str = region.get('segment')

    segment_num, gene_name = extract_segment_number_and_gene(segment_str)

    return segment_num if segment_num is not None else 999

def process_influenza_files(temp_folder, influenza_type, min_depth, min_cov):
    """
    Process Genome Detective assignment JSON files for influenza samples.
    Handles both old and new JSON formats.
    """
    segment_table_rows = []
    gisaid_table_rows = []
    
    # Track processed samples using the unique ID
    processed_samples = set()

    for json_file in Path(temp_folder).glob('*.assignments.json'):
        influenza_data_list = parse_json_file(json_file, influenza_type, 'influenza', min_depth, min_cov)
        if not influenza_data_list:
            continue

        for influenza_data in influenza_data_list:
            # Create unique key for this sample using the normalized ID and virus type
            sample_key = f"{influenza_data['LimsID']}_{influenza_data['virusType']}"
            
            # Skip if we've already processed this sample
            if sample_key in processed_samples:
                logging.debug(f"Skipping duplicate sample: {sample_key}")
                continue
                
            processed_samples.add(sample_key)
            
            # Skip if virus type is unknown
            if influenza_data['virusType'] not in ['A', 'B']:
                continue
                
            # Sort regions numerically
            sorted_regions = sorted(
                influenza_data['regions'], 
                key=sort_segment_key
            )

            # Collect segment metrics
            for region in sorted_regions:
                segment_num = region.get('segment_num')
                if segment_num is None:
                    continue
                    
                gene_name = VIRUS_CONFIG['influenza']['types'][
                    influenza_data['virusType']
                ]['segments'].get(segment_num, f"Segment{segment_num}")
                
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

            # Check if sample is a control
            is_control = "Yes" if any(ctrl in influenza_data['LimsID'] for ctrl in ["PC", "NC", "Neg", "Pos", "ERCC"]) else "No"

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
    Process Genome Detective assignment JSON files for non-segmented viruses.
    Handles both old and new JSON formats.
    """
    gisaid_table_rows = []
    
    # Track processed samples using the unique ID
    processed_samples = set()

    for json_file in Path(temp_folder).glob('*.assignments.json'):
        virus_data_list = parse_json_file(json_file, '', virus_key, min_depth, min_cov)
        if not virus_data_list:
            continue

        for virus_data in virus_data_list:
            # Skip duplicates using the normalized ID
            if virus_data['LimsID'] in processed_samples:
                logging.debug(f"Skipping duplicate sample: {virus_data['LimsID']}")
                continue
                
            processed_samples.add(virus_data['LimsID'])
            
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
    Handles both old and new JSON formats.
    """
    full_genome_data = []
    
    # Track processed samples using the unique ID
    processed_samples = set()

    for json_file in Path(temp_folder).glob('*.assignments.json'):
        try:
            with open(json_file, 'r') as f:
                data = json.load(f)

            base_name = os.path.basename(json_file).replace('.assignments.json', '')
            lims_id = normalize_limsid(base_name)
            
            strains = data.get('data', {}).get('attributes', {}).get('strains', [])

            for strain in strains:
                taxonomy = strain.get('taxonomyName', '').lower()
                
                # Create unique key using LIMS ID and taxonomy
                sample_key = f"{lims_id}_{taxonomy}"
                
                if sample_key in processed_samples:
                    continue
                    
                processed_samples.add(sample_key)

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
                    matched = False
                    for kw in VIRUS_CONFIG[virus_key]['taxonomy_keywords']:
                        if kw.lower() in taxonomy:
                            matched = True
                            break
                    
                    if matched:
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
    """
    print("\nAre the samples from multiple labs?")
    choice = input("Enter [y/n]: ").lower().strip()
    while choice not in ['y', 'n']:
        print("Invalid choice. Please enter y or n.")
        choice = input("Enter [y/n]: ").lower().strip()
    return choice == 'y'

def get_lab_info(multiple_labs, all_lims_ids):
    """
    Collect lab assignment information for each sample.
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
                parts = entry.split()
                if len(parts) < 2:
                    print("Invalid format. Use: LabName StartID-EndID")
                    continue

                lab_name = parts[0]
                range_str = parts[1]
                start, end = range_str.split('-')

                # Pattern to match C0XXXX or K0XXXX
                pattern = re.compile(r'([CK]0)(\d+)')
                start_match = pattern.match(start)
                end_match = pattern.match(end)

                if not start_match or not end_match:
                    print("Invalid ID format. Use format like C0001-C0010 or K001-K010")
                    continue

                prefix = start_match.group(1)  # This will be "C0" or "K0"
                start_num = int(start_match.group(2))
                end_num = int(end_match.group(2))

                matched_ids = []
                for lims_id in remaining:
                    # Extract the numeric part from the LIMS ID
                    lims_match = pattern.search(lims_id)
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
    """
    if lab_info is None:
        multiple_labs = prompt_for_multiple_labs()
        lab_info = get_lab_info(multiple_labs, all_lims_ids)

    if virus_key == 'influenza':
        if virus_type == 'C':
            default_header = "><isolate>/South Africa/<lab>-CERI-<LimsID>/2026_<gene>"
        else:
            default_header = f">{virus_type}/South Africa/<lab>-CERI-<LimsID>/2026_<gene>"
    elif virus_key == 'rsv':
        default_header = f">h{VIRUS_CONFIG[virus_key]['common_name']}/A/South Africa/<lab>-CERI-<LimsID>/2026"
    else:
        default_header = f">{VIRUS_CONFIG[virus_key]['common_name']}/South Africa/<lab>-CERI-<LimsID>/2026"

    print(f"\nDefault header format: {default_header}")
    confirm = input("Are you happy with this format? (y/n): ").lower()

    if confirm != 'y':
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
    """
    if len(lab_info) == 1:
        return next(iter(lab_info.keys()))

    for lab_name, lims_ids in lab_info.items():
        if limsid in lims_ids:
            return lab_name

    return "UNKNOWN"

def get_segment_number_and_gene(header_line, virus_name, virus_type):
    """
    Extract segment number and corresponding gene name from a segment label.
    """
    try:
        match = re.search(r'[Ss]egment\s*(\d+)', str(header_line))
        if match:
            segment_num = int(match.group(1))
            gene_name = VIRUS_CONFIG[virus_name]['types'][virus_type]['segments'].get(
                segment_num, f"Segment{segment_num}"
            )
            return segment_num, gene_name

        elif str(header_line).strip().isdigit():
            segment_num = int(header_line.strip())
            gene_name = VIRUS_CONFIG[virus_name]['types'][virus_type]['segments'].get(
                segment_num, f"Segment{segment_num}"
            )
            return segment_num, gene_name

    except Exception as e:
        logging.warning(f"Error parsing segment from header '{header_line}': {e}")

    return None, None

def find_sample_folder(input_folder, limsid):
    """
    Find the sample folder that contains the given LIMS ID.
    Uses the unique C0/K0 ID to match folders that may have prefixes/suffixes.
    """
    unique_id = extract_unique_id(limsid)
    
    for folder in os.listdir(input_folder):
        folder_path = os.path.join(input_folder, folder)
        if not os.path.isdir(folder_path):
            continue
            
        # Check if the unique ID is in the folder name
        if unique_id in folder:
            return folder_path
        
        # Also check if the full LIMS ID is in the folder name
        if limsid in folder:
            return folder_path
    
    return None

def find_alignment_file(sample_folder, ref_id):
    """
    Find alignment file in sample folder that matches the reference ID.
    Handles various filename patterns.
    """
    if not sample_folder or not os.path.exists(sample_folder):
        return None
        
    for file in os.listdir(sample_folder):
        file_path = os.path.join(sample_folder, file)
        if not os.path.isfile(file_path):
            continue
            
        # Check if reference ID is in the filename and it's an alignment file
        if ref_id in file and file.endswith('alignment-nt.fasta'):
            return file_path
        
        # Also try with just the accession number (without version)
        accession = ref_id.split('.')[0] if '.' in ref_id else ref_id
        if accession in file and file.endswith('alignment-nt.fasta'):
            return file_path
    
    return None

def write_segment_fastas(input_folder, segment_df, gisaid_df, virus_type):
    """
    Write individual FASTA files for each segment of influenza virus samples that passed GISAID thresholds.
    """
    eligible = gisaid_df[gisaid_df['Submit to GISAID'] == 'Yes']
    if eligible.empty:
        logging.info("No eligible samples for FASTA generation")
        return
        
    all_lims_ids = eligible['LimsID'].tolist()

    temp2_folder = os.path.join(os.getcwd(), 'Outputs')
    os.makedirs(temp2_folder, exist_ok=True)

    header_template, lab_info = prompt_user_for_header_format(virus_type, 'influenza', all_lims_ids)
    
    # Track processed samples to avoid duplicates
    processed_samples = set()

    for _, row in eligible.iterrows():
        limsid = row['LimsID']
        unique_id = extract_unique_id(limsid)
        
        # Skip duplicates
        sample_key = f"{limsid}_{row['VirusType']}"
        if sample_key in processed_samples:
            logging.debug(f"Skipping duplicate sample: {sample_key}")
            continue
            
        processed_samples.add(sample_key)

        # Determine output folder - use the unique ID for folder name
        if virus_type == 'C':
            dst_folder = os.path.join(temp2_folder, unique_id, f'Influenza_{row["VirusType"]}')
        else:
            dst_folder = os.path.join(temp2_folder, unique_id)
        os.makedirs(dst_folder, exist_ok=True)

        # Find sample input folder
        sample_input_folder = find_sample_folder(input_folder, limsid)

        if not sample_input_folder:
            logging.warning(f"No folder containing {limsid} (unique ID: {unique_id}) found in {input_folder}")
            continue

        # Get segments for this sample
        segments = segment_df[(segment_df['LimsID'] == limsid) & 
                              (segment_df['VirusType'] == row['VirusType'])]

        for _, seg_row in segments.iterrows():
            ref_id = seg_row['ReferenceSequenceId']
            depth_id = seg_row['DepthOfCoverage']
            coverage_id = seg_row['CoveragePercentage']
            seg_name = seg_row['Segment']
            seg_num, gene_name = get_segment_number_and_gene(seg_name, 'influenza', row['VirusType'])

            if not seg_num:
                logging.warning(f"Could not extract segment number from {seg_name} for sample {limsid}")
                continue

            if depth_id >= 10 and coverage_id >= 80:
                alignment_file = find_alignment_file(sample_input_folder, ref_id)
                
                if alignment_file:
                    try:
                        with open(alignment_file, 'r') as f:
                            lines = [line.strip() for line in f if line.strip()]

                        headers = [i for i, line in enumerate(lines) if line.startswith('>')]
                        if len(headers) != 2:
                            logging.warning(f"Expected exactly 2 sequences in {alignment_file}, found {len(headers)}. Skipping.")
                            continue

                        consensus_start = headers[1]
                        full_sequence = ''.join(lines[consensus_start + 1:])

                        lab_name = get_lab_for_limsid(limsid, lab_info)

                        new_header = header_template \
                            .replace('<LimsID>', unique_id) \
                            .replace('<gene>', gene_name) \
                            .replace('<lab>', lab_name) \
                            .replace('<isolate>', f'{row["VirusType"]}')

                        if not new_header.startswith('>'):
                            new_header = '>' + new_header

                        new_header = new_header.replace(f'>{row["VirusType"]}/{row["VirusType"]}/', f'>{row["VirusType"]}/')
                        new_header = new_header.replace('>AB/', f'>{row["VirusType"]}/')

                        dst_file = os.path.join(dst_folder, f'segment{seg_num}.fasta')
                        
                        # Check if file already exists
                        if os.path.exists(dst_file):
                            logging.debug(f"File already exists, skipping: {dst_file}")
                            continue
                            
                        with open(dst_file, 'w') as out:
                            out.write(f"{new_header}\n{full_sequence}\n")

                        logging.info(f"Written segment fasta: {dst_file}")

                    except Exception as e:
                        logging.error(f"Failed to process {alignment_file}: {e}")
                else:
                    logging.debug(f"No alignment file found for reference ID {ref_id} in {sample_input_folder}")

def write_virus_fastas(input_folder, df, virus_key, lab_info):
    """
    Write FASTA files for non-segmented virus samples that passed GISAID criteria.
    """
    temp_folder = os.path.join(os.getcwd(), 'Outputs')
    os.makedirs(temp_folder, exist_ok=True)

    print(f"\nDefault header format for {VIRUS_CONFIG[virus_key]['name']}:")
    if virus_key == "rsv":
        default_header = f">h{VIRUS_CONFIG[virus_key]['common_name']}/A/South Africa/<lab>-CERI-<LimsID>/2026"
    else:
         default_header = f">{VIRUS_CONFIG[virus_key]['common_name']}/South Africa/<lab>-CERI-<LimsID>/2026"   
    print(default_header)

    confirm = input("Use this format? (y/n): ").lower()
    if confirm != 'y':
        header_template = input("Enter custom format using <LimsID> and <lab>: ")
    else:
        header_template = default_header

    eligible_samples = df[df['Submit to GISAID'] == 'Yes']
    if eligible_samples.empty:
        logging.info("No eligible samples for FASTA generation")
        return
    
    # Track processed samples
    processed_samples = set()

    for _, row in eligible_samples.iterrows():
        limsid = row['LimsID']
        unique_id = extract_unique_id(limsid)
        
        # Skip duplicates
        if limsid in processed_samples:
            logging.debug(f"Skipping duplicate sample: {limsid}")
            continue
            
        processed_samples.add(limsid)
        
        ref_id = row['ReferenceSequenceId']

        sample_folder = os.path.join(temp_folder, unique_id)
        os.makedirs(sample_folder, exist_ok=True)

        # Find sample input folder
        sample_input_folder = find_sample_folder(input_folder, limsid)

        if not sample_input_folder:
            logging.warning(f"No folder containing {limsid} (unique ID: {unique_id}) found in {input_folder}")
            continue

        alignment_file = find_alignment_file(sample_input_folder, ref_id)
        
        if alignment_file:
            try:
                with open(alignment_file, 'r') as f:
                    lines = [line.strip() for line in f if line.strip()]

                headers = [i for i, line in enumerate(lines) if line.startswith('>')]
                if len(headers) < 2:
                    logging.warning(f"Expected 2 headers in {alignment_file}, found {len(headers)}. Skipping.")
                    continue

                sequence = ''.join(lines[headers[1] + 1:])

                lab_name = get_lab_for_limsid(limsid, lab_info)
                new_header = header_template\
                    .replace('<LimsID>', unique_id)\
                    .replace('<lab>', lab_name)

                dst_file = os.path.join(sample_folder, f"{unique_id}_{virus_key}.fasta")
                
                # Check if file already exists
                if os.path.exists(dst_file):
                    logging.debug(f"File already exists, skipping: {dst_file}")
                    continue
                    
                with open(dst_file, 'w') as out:
                    out.write(f"{new_header}\n{sequence}\n")

                logging.info(f"Written FASTA file: {dst_file}")

            except Exception as e:
                logging.error(f"Failed to process {alignment_file}: {e}")
        else:
            logging.debug(f"No alignment file found for reference ID {ref_id} in {sample_input_folder}")

def concatenate_segments_per_sample(virus_type):
    """
    Concatenate individual influenza segment FASTA files per sample into one multi-segment file.
    """
    temp2_folder = os.path.join(os.getcwd(), 'Outputs')
    
    # Track processed samples
    processed_samples = set()

    for folder_name in os.listdir(temp2_folder):
        sample_folder = os.path.join(temp2_folder, folder_name)
        if not os.path.isdir(sample_folder):
            continue
        
        # The folder name should be the unique ID
        unique_id = folder_name

        if virus_type == 'AB':
            for influ_type in ['A', 'B']:
                type_folder = os.path.join(sample_folder, f'Influenza_{influ_type}')
                if not os.path.exists(type_folder):
                    continue
                    
                # Create unique key
                sample_key = f"{unique_id}_{influ_type}"
                if sample_key in processed_samples:
                    continue
                    
                processed_samples.add(sample_key)

                segment_files = sorted([
                    os.path.join(type_folder, f)
                    for f in os.listdir(type_folder)
                    if f.startswith('segment') and f.endswith('.fasta')
                ])

                if segment_files:
                    concatenated_fasta = os.path.join(type_folder, f"{unique_id}_all_segments.fasta")
                    
                    # Skip if already exists
                    if os.path.exists(concatenated_fasta):
                        logging.debug(f"Concatenated file already exists: {concatenated_fasta}")
                        continue
                        
                    with open(concatenated_fasta, 'w') as outfile:
                        for seg_file in segment_files:
                            with open(seg_file, 'r') as infile:
                                content = infile.read().strip()
                                outfile.write(content + '\n')
                    #logging.info(f"Concatenated segments for sample {unique_id} (type {influ_type}) into {concatenated_fasta}")

        else:
            # Single-type influenza sample
            if unique_id in processed_samples:
                continue
                
            processed_samples.add(unique_id)
            
            segment_files = sorted([
                os.path.join(sample_folder, f)
                for f in os.listdir(sample_folder)
                if f.startswith('segment') and f.endswith('.fasta')
            ])

            if segment_files:
                concatenated_fasta = os.path.join(sample_folder, f"{unique_id}_all_segments.fasta")
                
                # Skip if already exists
                if os.path.exists(concatenated_fasta):
                    logging.debug(f"Concatenated file already exists: {concatenated_fasta}")
                    continue
                    
                with open(concatenated_fasta, 'w') as outfile:
                    for seg_file in segment_files:
                        with open(seg_file, 'r') as infile:
                            content = infile.read().strip()
                            outfile.write(content + '\n')
                logging.info(f"Concatenated segments for sample {unique_id} into {concatenated_fasta}")

def concatenate_all_samples_fasta(virus_type, virus_key, output_filename_prefix="all_samples_combined"):
    """
    Concatenate all per-sample full-genome FASTA files into a single master file.
    """
    temp2_folder = os.path.join(os.getcwd(), 'Outputs')

    if virus_type == 'AB':
        for influ_type in ['A', 'B']:
            combined_fasta_path = os.path.join(temp2_folder, f"{output_filename_prefix}_{influ_type}.fasta")
            
            # Skip if already exists
            if os.path.exists(combined_fasta_path):
                logging.info(f"Combined file already exists: {combined_fasta_path}")
                continue
                
            with open(combined_fasta_path, 'w') as outfile:
                for folder_name in os.listdir(temp2_folder):
                    sample_folder = os.path.join(temp2_folder, folder_name)
                    if not os.path.isdir(sample_folder):
                        continue

                    if any(ctrl in folder_name for ctrl in ["PC", "NC", "Neg", "Pos", "ERCC"]):
                        logging.info(f"Skipping control sample: {folder_name}")
                        continue

                    type_folder = os.path.join(sample_folder, f'Influenza_{influ_type}')
                    if not os.path.exists(type_folder):
                        continue

                    all_segments_fasta = os.path.join(type_folder, f"{folder_name}_all_segments.fasta")
                    if not os.path.exists(all_segments_fasta):
                        logging.debug(f"Missing concatenated fasta for sample {folder_name}: {all_segments_fasta}")
                        continue

                    with open(all_segments_fasta, 'r') as infile:
                        content = infile.read().strip()
                        outfile.write(f">{folder_name}\n")
                        outfile.write(content + '\n')

            logging.info(f"Combined fasta for all {influ_type} samples written to {combined_fasta_path}")

    else:
        combined_fasta_path = os.path.join(temp2_folder, f"{output_filename_prefix}.fasta")
        
        # Skip if already exists
        if os.path.exists(combined_fasta_path):
            logging.info(f"Combined file already exists: {combined_fasta_path}")
            return
            
        with open(combined_fasta_path, 'w') as outfile:
            for folder_name in os.listdir(temp2_folder):
                sample_folder = os.path.join(temp2_folder, folder_name)
                if not os.path.isdir(sample_folder):
                    continue

                if any(ctrl in folder_name for ctrl in ["PC", "NC", "Neg", "Pos", "ERCC"]):
                    logging.info(f"Skipping control sample: {folder_name}")
                    continue

                # Look for the concatenated FASTA file
                if virus_key == 'influenza':
                    fasta_file = os.path.join(sample_folder, f"{folder_name}_all_segments.fasta")
                else:
                    # For non-influenza viruses, look for the virus-specific file
                    fasta_file = os.path.join(sample_folder, f"{folder_name}_{virus_key}.fasta")
                    
                if not os.path.exists(fasta_file):
                    logging.debug(f"Missing fasta file for sample {folder_name}: {fasta_file}")
                    continue

                with open(fasta_file, 'r') as infile:
                    content = infile.read().strip()
                    outfile.write(content + '\n')

        logging.info(f"Combined fasta for all samples written to {combined_fasta_path}")

def main(input_folder):
    """
    Main entry point for the viral analysis pipeline.
    """
    logging.info("Starting analysis")

    virus_type, virus_key = get_virus_type()
    virus_type_header = (
        'A' if virus_type == 'A'
        else 'B' if virus_type == 'B'
        else 'AB' if virus_key == 'influenza'
        else virus_key
    )
    min_depth, min_cov = get_coverage_thresholds()

    print("Copying assignment files to temporary folder...")
    temp_folder = create_temp_folder(input_folder)

    if virus_key == 'influenza':
        segment_df, gisaid_df = process_influenza_files(temp_folder, virus_type, min_depth, min_cov)
    else:
        gisaid_df = process_non_segmented_virus(temp_folder, virus_key, min_depth, min_cov)
        segment_df = None

    full_genome_df = extract_genome_info(temp_folder, virus_type, virus_key)

    # Validation checks
    if virus_key == 'influenza':
        if virus_type in ['A', 'B'] and (gisaid_df.empty or len(gisaid_df) == 0):
            logging.error(f"Your chosen virus type (Influenza {virus_type}) was not detected in any samples")
            exit(1)
        elif virus_type == 'C' and not gisaid_df.empty and (
            len(gisaid_df[gisaid_df['VirusType'] == 'A']) == 0 and
            len(gisaid_df[gisaid_df['VirusType'] == 'B']) == 0
        ):
            logging.error("Neither Influenza A nor B were detected in any samples")
            exit(1)
    else:
        if gisaid_df.empty or len(gisaid_df) == 0:
            logging.error(f"Your chosen virus ({VIRUS_CONFIG[virus_key]['name']}) was not detected in any samples")
            exit(1)

    # Create output directory for summaries
    output_dir = os.path.join(os.getcwd(), 'Summary_files')
    os.makedirs(output_dir, exist_ok=True)

    # Generate summary CSVs
    if virus_key == 'influenza':
        if virus_type == 'C':
            for v_type in ['A', 'B']:
                type_segment_df = segment_df[segment_df['VirusType'] == v_type].sort_values(by='LimsID') if not segment_df.empty else pd.DataFrame()
                type_gisaid_df = gisaid_df[gisaid_df['VirusType'] == v_type].sort_values(by='LimsID') if not gisaid_df.empty else pd.DataFrame()
                type_genome_df = full_genome_df[full_genome_df['VirusType'] == v_type].sort_values(by='LimsID') if not full_genome_df.empty else pd.DataFrame()

                if not type_segment_df.empty:
                    type_segment_df.to_csv(os.path.join(output_dir, f'influenza_{v_type.lower()}_segments.csv'), index=False)
                if not type_gisaid_df.empty:
                    type_gisaid_df.to_csv(os.path.join(output_dir, f'gisaid_submission_status_{v_type.lower()}.csv'), index=False)
                if not type_genome_df.empty:
                    type_genome_df.to_csv(os.path.join(output_dir, f'full_genome_info_{v_type.lower()}.csv'), index=False)

                if not type_gisaid_df.empty:
                    eligible_count = len(type_gisaid_df[type_gisaid_df['Submit to GISAID'] == 'Yes'])
                    control_count = len(type_gisaid_df[(type_gisaid_df['Submit to GISAID'] == 'Yes') &
                                                       (type_gisaid_df['Is Control'] == 'Yes')])
                    logging.info(f"Influenza {v_type} - Total samples: {len(type_gisaid_df)}, "
                                 f"GISAID-eligible: {eligible_count} (including {control_count} controls)")
        else:
            if not segment_df.empty:
                segment_df = segment_df.sort_values(by='LimsID')
                segment_df.to_csv(os.path.join(output_dir, f'influenza_{virus_type.lower()}_segments.csv'), index=False)
            
            if not gisaid_df.empty:
                gisaid_df = gisaid_df.sort_values(by='LimsID')
                gisaid_df.to_csv(os.path.join(output_dir, f'gisaid_submission_status_{virus_type.lower()}.csv'), index=False)
            
            if not full_genome_df.empty:
                full_genome_df = full_genome_df.sort_values(by='LimsID')
                full_genome_df.to_csv(os.path.join(output_dir, f'full_genome_info_{virus_type.lower()}.csv'), index=False)

            if not gisaid_df.empty:
                eligible_count = len(gisaid_df[gisaid_df['Submit to GISAID'] == 'Yes'])
                control_count = len(gisaid_df[(gisaid_df['Submit to GISAID'] == 'Yes') &
                                              (gisaid_df['Is Control'] == 'Yes')])
                logging.info(f"Total samples: {len(gisaid_df)}, GISAID-eligible: {eligible_count} (including {control_count} controls)")
    else:
        if not gisaid_df.empty:
            gisaid_df = gisaid_df.sort_values(by='LimsID')
            gisaid_df.to_csv(os.path.join(output_dir, f'gisaid_submission_status_{virus_key}.csv'), index=False)
        
        if not full_genome_df.empty:
            full_genome_df = full_genome_df.sort_values(by='LimsID')
            full_genome_df.to_csv(os.path.join(output_dir, f'full_genome_info_{virus_key}.csv'), index=False)

        if not gisaid_df.empty:
            eligible_count = len(gisaid_df[gisaid_df['Submit to GISAID'] == 'Yes'])
            control_count = len(gisaid_df[(gisaid_df['Submit to GISAID'] == 'Yes') &
                                          (gisaid_df['Is Control'] == 'Yes')])
            logging.info(f"Total samples: {len(gisaid_df)}, GISAID-eligible: {eligible_count} (including {control_count} controls)")

    # Write FASTA files
    if virus_key == 'influenza' and not gisaid_df.empty:
        logging.info("Writing segment FASTA files...")
        write_segment_fastas(input_folder, segment_df, gisaid_df, virus_type)
        logging.info("Concatenating segments per sample...")
        concatenate_segments_per_sample(virus_type_header)
        logging.info("Creating combined FASTA file...")
        concatenate_all_samples_fasta(virus_type_header, 'influenza')
    elif not gisaid_df.empty:
        logging.info("Preparing lab information...")
        all_lims_ids = gisaid_df['LimsID'].tolist()
        multiple_labs = prompt_for_multiple_labs()
        lab_info = get_lab_info(multiple_labs, all_lims_ids)
        logging.info("Writing virus FASTA files...")
        write_virus_fastas(input_folder, gisaid_df, virus_key, lab_info)
        logging.info("Creating combined FASTA file...")
        concatenate_all_samples_fasta(virus_type_header, virus_key)

    logging.info("Pipeline completed successfully!")
    return True

if __name__ == "__main__":
    import argparse
    import logging

    parser = argparse.ArgumentParser(
        description=' Viral GISAID Submission Pipeline\n'
                    'Processes Genome Detective outputs and prepares submission files.'
    )
    parser.add_argument(
        'input_folder',
        help='Path to input folder containing LimsID-named subdirectories with .assignments.json and alignment FASTAs'
    )
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s [%(levelname)s] %(message)s',
        handlers=[logging.StreamHandler()]
    )

    success = main(args.input_folder)

    if not success:
        logging.error("Pipeline failed.")
        exit(1)