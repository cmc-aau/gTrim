use assert_cmd::cargo::cargo_bin_cmd;
use std::path::Path;
use std::process::Output;
use tempfile::tempdir;

mod test_utils;
use test_utils::{compare_fasta_files, compare_mask_files, compare_tsv_files};

const MARKER_INFO: &str = "./tests/input_file/trim_test_markers_info.tsv";

fn run_trim(alignment: &str, extra_args: &[&str], output_dir: &Path) -> Output {
    let mut cmd = cargo_bin_cmd!("gtrim");
    cmd.arg("trim")
        .arg("-i")
        .arg(alignment)
        .arg("-m")
        .arg(MARKER_INFO)
        .args(extra_args)
        .arg("-o")
        .arg(output_dir);
    cmd.output().expect("failed to execute process")
}

fn assert_success(output: &Output) {
    if !output.status.success() {
        println!("stdout: {}", String::from_utf8_lossy(&output.stdout));
        println!("stderr: {}", String::from_utf8_lossy(&output.stderr));
        panic!("Command executed with failing error code");
    }
}

fn compare_outputs(output_dir: &Path, prefix: &str) {
    let generated = |suffix: &str| output_dir.join(format!("{}{}", prefix, suffix));
    let expected = |suffix: &str| Path::new("./tests/ground_truth").join(format!("{}{}", prefix, suffix));

    compare_mask_files(&generated(".filtered.mask"), &expected(".filtered.mask"));
    compare_mask_files(&generated(".mask"), &expected(".mask"));
    compare_tsv_files(&generated(".markers.tsv"), &expected(".markers.tsv"), 1e-6, 10);
    compare_tsv_files(&generated(".genomes.tsv"), &expected(".genomes.tsv"), 1e-6, 4);
    compare_fasta_files(&generated(".fasta.gz"), &expected(".fasta.gz"));
}

// Test data: the first 200 genomes and first 10 bac120 markers (1,792 columns) of
// test_sets/gtdb_n5000_test/gtdb_r226_n5000.faa.gz. The filtered mask and marker statistics
// match GenomeTreeTk `msa_filter` (no WitChi file) with --msa_length 500. Sampled columns use the
// default seed (42).
#[test]
fn test_trim() {
    let temp_dir = tempdir().unwrap();
    let output = run_trim("./tests/input_file/trim_test.faa.gz", &["--msa-length", "500", "-t", "2"], temp_dir.path());
    assert_success(&output);
    compare_outputs(temp_dir.path(), "trim_test.trimmed");
}

// trim_test_witchi.faa.gz is trim_test.faa.gz with the 150 columns listed in trim_test_witchi.tsv
// (a WitChi pruned TSV) removed; trim_test_witchi.mask is the same mask in binary form. None of the
// removed columns is the first column of a marker. The filtered mask, marker statistics and
// number of columns sampled from each marker match GenomeTreeTk `msa_filter` given the WitChi file,
// with --msa_length 500.
#[test]
fn test_trim_witchi_mask() {
    let temp_dir = tempdir().unwrap();
    let output = run_trim(
        "./tests/input_file/trim_test_witchi.faa.gz",
        &["--mask-file", "./tests/input_file/trim_test_witchi.tsv", "--msa-length", "500", "-t", "2"],
        temp_dir.path(),
    );
    assert_success(&output);
    compare_outputs(temp_dir.path(), "trim_test_witchi.trimmed");
}

#[test]
fn test_trim_binary_mask() {
    let temp_dir = tempdir().unwrap();
    let output = run_trim(
        "./tests/input_file/trim_test_witchi.faa.gz",
        &["--mask-file", "./tests/input_file/trim_test_witchi.mask", "--msa-length", "500"],
        temp_dir.path(),
    );
    assert_success(&output);
    compare_outputs(temp_dir.path(), "trim_test_witchi.trimmed");
}

#[test]
fn test_trim_length_mismatch() {
    let temp_dir = tempdir().unwrap();
    let output = run_trim("./tests/input_file/gap_test.fna.gz", &[], temp_dir.path());
    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("does not equal the summed length"));
}

#[test]
fn test_trim_mask_mismatch() {
    // unmasked alignment given with a mask that removes columns
    let temp_dir = tempdir().unwrap();
    let output = run_trim(
        "./tests/input_file/trim_test.faa.gz",
        &["--mask-file", "./tests/input_file/trim_test_witchi.tsv"],
        temp_dir.path(),
    );
    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("does not equal the number of columns retained by the mask file"));

    // mask that does not match the marker info file
    let output = run_trim(
        "./tests/input_file/trim_test_witchi.faa.gz",
        &["--mask-file", "./tests/ground_truth/gap_test.pruned.mask"],
        temp_dir.path(),
    );
    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("summed length of 1792 columns"));
}
