use std::fs;
use std::path::Path;
use std::collections::HashMap;

use anyhow::{Context, Result};
use needletail::parse_fastx_file;

pub fn compare_tsv_files(generated_path: &Path, ground_truth_path: &Path, tolerance: f64, cols_to_compare: usize) {
    let generated_content = fs::read_to_string(generated_path)
        .unwrap_or_else(|err| panic!("Could not read generated file {:?}: {}", generated_path, err));
    let ground_truth_content = fs::read_to_string(ground_truth_path).expect("Could not read ground truth file");

    let generated_lines: Vec<&str> = generated_content.lines().collect();
    let ground_truth_lines: Vec<&str> = ground_truth_content.lines().collect();

    assert_eq!(generated_lines.len(), ground_truth_lines.len(), "Line count differs between {:?} and {:?}", generated_path, ground_truth_path);

    for (i, (gen_line, truth_line)) in generated_lines.iter().zip(ground_truth_lines.iter()).enumerate() {
        let gen_fields: Vec<&str> = gen_line.split('\t').collect();
        let truth_fields: Vec<&str> = truth_line.split('\t').collect();

        for j in 0..cols_to_compare {
            let gen_field = gen_fields[j];
            let truth_field = truth_fields[j];
            if let (Ok(gen_float), Ok(truth_float)) = (gen_field.parse::<f64>(), truth_field.parse::<f64>()) {
                assert!((gen_float - truth_float).abs() < tolerance, "Floating point difference on line {}, field {} between {:?} and {:?}. Gen: {}, Truth: {}", i + 1, j + 1, generated_path, ground_truth_path, gen_float, truth_float);
            } else {
                assert_eq!(gen_field, truth_field, "Field difference on line {}, field {} between {:?} and {:?}. Gen: '{}', Truth: '{}'", i + 1, j + 1, generated_path, ground_truth_path, gen_field, truth_field);
            }
        }
    }
}

fn read_fasta_to_map(path: &Path) -> Result<HashMap<String, String>> {
    let mut reader = parse_fastx_file(path)
            .with_context(|| format!("Failed to open alignment file: {}", path.display()))?;

    let mut fasta_map = HashMap::new();

    while let Some(record) = reader.next() {
        let rec = record.with_context(|| "Failed to parse FASTA record")?;
        let id = String::from_utf8_lossy(rec.id()).to_string();
        let seq: String = rec.seq().iter().map(|&b| (b as char).to_ascii_uppercase()).collect();

        fasta_map.insert(id, seq);
    }

    Ok(fasta_map)
}

pub fn compare_fasta_files(generated_path: &Path, ground_truth_path: &Path)  {
    let generated_map = read_fasta_to_map(generated_path).unwrap_or_else(|err| panic!("Could not read FASTA file {:?}: {}", generated_path, err));
    let ground_truth_map = read_fasta_to_map(ground_truth_path).unwrap_or_else(|err| panic!("Could not read FASTA file {:?}: {}", generated_path, err));

    assert_eq!(generated_map.len(), ground_truth_map.len(), "FASTA record count differs between {:?} and {:?}", generated_path, ground_truth_path);

    for (id, gen_seq) in generated_map.iter() {
        let truth_seq = ground_truth_map.get(id).unwrap_or_else(|| panic!("Sequence ID '{}' not found in ground truth file {:?}", id, ground_truth_path));
        assert_eq!(gen_seq, truth_seq, "Sequence mismatch for ID '{}' between {:?} and {:?}", id, generated_path, ground_truth_path);
    }
}

#[allow(dead_code)]
pub fn compare_mask_files(generated_path: &Path, ground_truth_path: &Path)  {
    let generated_mask = fs::read_to_string(generated_path).unwrap_or_else(|err| panic!("Could not read generated mask {:?}: {}", generated_path, err));
    let ground_truth_mask = fs::read_to_string(ground_truth_path).expect("Could not read ground truth mask");

    let generated_mask: Vec<&str> = generated_mask.lines().collect();
    let ground_truth_mask: Vec<&str> = ground_truth_mask.lines().collect();

    // check that masks are identical
    assert_eq!(generated_mask.len(), ground_truth_mask.len(), "Mask line count differs between {:?} and {:?}", generated_path, ground_truth_path);
    for (i, (gen_line, truth_line)) in generated_mask.iter().zip(ground_truth_mask.iter()).enumerate() {
        assert_eq!(gen_line, truth_line, "Mask content differs on line {} between {:?} and {:?}. Gen: '{}', Truth: '{}'", i + 1, generated_path, ground_truth_path, gen_line, truth_line);
    }
}