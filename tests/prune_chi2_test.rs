use assert_cmd::cargo::cargo_bin_cmd;
use std::path::Path;
use tempfile::tempdir;

mod test_utils;
use test_utils::{compare_tsv_files, compare_fasta_files};

#[test]
fn test_prune_chi2() {
    let temp_dir = tempdir().unwrap();
    let output_dir = temp_dir.path();

    let mut cmd = cargo_bin_cmd!("gtrim");
    cmd.arg("prune-chi2")
        .arg("--file")
        .arg("./tests/input_file/gtdb_r226-n100.faa.gz")
        .arg("--residues-to-prune-perc")
        .arg("0")
        .arg("--max-residue-pruned")
        .arg("10")
        .arg("--permutations")
        .arg("10")
        .arg("--out-dir")
        .arg(output_dir);

    let output = cmd.output().expect("failed to execute process");
    if !output.status.success() {
        println!("stdout: {}", String::from_utf8_lossy(&output.stdout));
        println!("stderr: {}", String::from_utf8_lossy(&output.stderr));
        panic!("Command executed with failing error code");
    }

    let generated_stats_path = output_dir.join("gtdb_r226-n100_p10_s1_m10_pruned.tsv");
    let expected_stats_path = Path::new("./tests/ground_truth/gtdb_r226-n100_quartic_s1_pruned.tsv");
    compare_tsv_files(&generated_stats_path, expected_stats_path, 1e-6, 5);

    let generated_fasta_path = output_dir.join("gtdb_r226-n100_p10_s1_m10_pruned.fasta.gz");
    let expected_fasta_path = Path::new("./tests/ground_truth/gtdb_r226-n100_quartic_s1_pruned.fasta");
    compare_fasta_files(&generated_fasta_path, expected_fasta_path);

    // Hard to compare these files as the output is highly dependent on the generated permutations.
    //let generated_scores_path = output_dir.join("gtdb_r226-n100_p10_s1_m10_pruned_scores.tsv");
    //let expected_scores_path = Path::new("./tests/ground_truth/gtdb_r226-n100_quartic_s1_pruned_scores.tsv");
    //compare_tsv_files(&generated_scores_path, expected_scores_path, 1e-6, 3);
}
