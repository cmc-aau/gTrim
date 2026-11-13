use assert_cmd::cargo::cargo_bin_cmd;
use std::path::Path;
use tempfile::tempdir;

mod test_utils;
use test_utils::{compare_fasta_files, compare_mask_files, compare_tsv_files};

// Test data: the first 200 genomes and first 10 bac120 markers (1,792 columns) of
// test_sets/gtdb_n5000_test/gtdb_r226_n5000.faa.gz. The filtered mask and marker statistics
// match GenomeTreeTk `msa_filter` (no WitChi file) with --msa_length 500. Sampled columns use the
// default seed (42).
#[test]
fn test_trim() {
    let temp_dir = tempdir().unwrap();
    let output_dir = temp_dir.path();

    let mut cmd = cargo_bin_cmd!("gtrim");
    cmd.arg("trim")
        .arg("-i")
        .arg("./tests/input_file/trim_test.faa.gz")
        .arg("-m")
        .arg("./tests/input_file/trim_test_markers_info.tsv")
        .arg("--msa-length")
        .arg("500")
        .arg("-t")
        .arg("2")
        .arg("-o")
        .arg(output_dir);

    let output = cmd.output().expect("failed to execute process");
    if !output.status.success() {
        println!("stdout: {}", String::from_utf8_lossy(&output.stdout));
        println!("stderr: {}", String::from_utf8_lossy(&output.stderr));
        panic!("Command executed with failing error code");
    }

    let prefix = "trim_test_l500_aa50_c25-95_trimmed";
    let generated = |suffix: &str| output_dir.join(format!("{}{}", prefix, suffix));
    let expected = |suffix: &str| Path::new("./tests/ground_truth").join(format!("{}{}", prefix, suffix));

    compare_mask_files(&generated("_filtered.mask"), &expected("_filtered.mask"));
    compare_mask_files(&generated(".mask"), &expected(".mask"));
    compare_tsv_files(&generated("_markers.tsv"), &expected("_markers.tsv"), 1e-6, 9);
    compare_tsv_files(&generated("_genomes.tsv"), &expected("_genomes.tsv"), 1e-6, 4);
    compare_fasta_files(&generated(".fasta.gz"), &expected(".fasta.gz"));
}

#[test]
fn test_trim_length_mismatch() {
    let temp_dir = tempdir().unwrap();

    let mut cmd = cargo_bin_cmd!("gtrim");
    cmd.arg("trim")
        .arg("-i")
        .arg("./tests/input_file/gap_test.fna.gz")
        .arg("-m")
        .arg("./tests/input_file/trim_test_markers_info.tsv")
        .arg("-o")
        .arg(temp_dir.path());

    let output = cmd.output().expect("failed to execute process");
    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("does not equal the summed length"));
}
