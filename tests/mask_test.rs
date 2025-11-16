use assert_cmd::cargo::cargo_bin_cmd;
use std::path::Path;
use tempfile::tempdir;

mod test_utils;
use test_utils::compare_fasta_files;


#[test]
fn test_prune_gap() {
    let temp_dir = tempdir().unwrap();
    let output_dir = temp_dir.path();

    let mut cmd = cargo_bin_cmd!("gtrim");
    cmd.arg("mask")
        .arg("-i")
        .arg("./tests/input_file/gap_test.fna.gz")
        .arg("--mask-file")
        .arg("./tests/ground_truth/gap_test_gt50_m100000_pruned.mask")
        .arg("-o")
        .arg(output_dir);

    let output = cmd.output().expect("failed to execute process");
    if !output.status.success() {
        println!("stdout: {}", String::from_utf8_lossy(&output.stdout));
        println!("stderr: {}", String::from_utf8_lossy(&output.stderr));
        panic!("Command executed with failing error code");
    }

    let generated_fasta_path = output_dir.join("gap_test_masked.fasta.gz");
    let expected_fasta_path = Path::new("./tests/ground_truth/gap_test_masked.fasta.gz");
    compare_fasta_files(&generated_fasta_path, expected_fasta_path);
}
