use assert_cmd::cargo::cargo_bin_cmd;
use std::fs;
use std::path::Path;
use tempfile::tempdir;

mod test_utils;
use test_utils::{compare_tsv_files, compare_fasta_files};

// Test data is the 6 shortest loci and 150 genomes from gtdb_r232_fun77_markers, created with:
//   tests/scripts/make_divergent_test_set.py <src_dir> <out_dir> --num-loci 6 --shortest --num-genomes 150 \
//       --seed 1 --plant 333279at4751:GCA_047655925.1:GCA_958449795.1
// The planted locus replaces the sequence of GCA_047655925.1 (~90% AAI to its closest neighbours)
// with the sequence of the distantly related GCA_958449795.1.
// The test set has only ~200-420 comparisons per locus, so a smaller running median window (40)
// than the default (120, intended for ~16k comparisons per locus) is used for both tests.
#[test]
fn test_mask_divergent_loci() {
    let temp_dir = tempdir().unwrap();
    let output_dir = temp_dir.path();

    let mut cmd = cargo_bin_cmd!("gtrim");
    cmd.arg("mask-divergent-loci")
        .arg("-i")
        .arg("./tests/input_file/divergent_test")
        .arg("-o")
        .arg(output_dir)
        .arg("--window-size")
        .arg("40")
        .arg("--upper-window-size")
        .arg("40")
        .arg("--taxonomy")
        .arg("./tests/input_file/divergent_test_taxonomy.tsv")
        .arg("--plot");

    let output = cmd.output().expect("failed to execute process");
    if !output.status.success() {
        println!("stdout: {}", String::from_utf8_lossy(&output.stdout));
        println!("stderr: {}", String::from_utf8_lossy(&output.stderr));
        panic!("Command executed with failing error code");
    }

    // planted divergent locus must be masked
    let masked = fs::read_to_string(output_dir.join("divergent_loci_masked.tsv")).unwrap();
    assert!(
        masked.lines().any(|line| line.starts_with("333279at4751\tGCA_047655925.1\t")),
        "Planted divergent locus was not masked:\n{}",
        masked
    );

    let truth_dir = Path::new("./tests/ground_truth/divergent_test");
    compare_tsv_files(
        &output_dir.join("divergent_loci_masked.tsv"),
        &truth_dir.join("divergent_loci_masked.tsv"),
        1e-6,
        13,
    );
    compare_tsv_files(
        &output_dir.join("divergent_loci_summary.tsv"),
        &truth_dir.join("divergent_loci_summary.tsv"),
        1e-6,
        11,
    );
    compare_tsv_files(
        &output_dir.join("divergent_loci_genome_summary.tsv"),
        &truth_dir.join("divergent_loci_genome_summary.tsv"),
        1e-6,
        11,
    );

    // planted locus is also identified by the high-identity test since the donor is in the data, its best
    // match is the donor, and the matching genomes form one lineage (synthetic taxonomy: same genus as the
    // donor, same phylum as the focal genome), while the donor genome itself is not masked
    let planted: Vec<&str> = masked
        .lines()
        .find(|line| line.starts_with("333279at4751\tGCA_047655925.1\t"))
        .unwrap()
        .split('\t')
        .collect();
    assert_eq!(planted[6], "both", "Planted divergent locus was not masked by both tests:\n{}", masked);
    assert_eq!(planted[7], "GCA_958449795.1", "Best match of planted locus is not the donor:\n{}", masked);
    assert_eq!(&planted[10..13], &["true", "g__Genus2", "p__Phylum1"], "Unexpected source of planted locus:\n{}", masked);
    assert!(!masked.lines().any(|line| line.split('\t').nth(1) == Some("GCA_958449795.1")), "Donor genome was masked:\n{}", masked);

    // the donor is the main source of the planted genome's high-identity masks
    let genome_summary = fs::read_to_string(output_dir.join("divergent_loci_genome_summary.tsv")).unwrap();
    assert!(
        genome_summary.lines().any(|line| line == "GCA_047655925.1\t6\t100.00\t1\t1\t1\t5\t83.33\tfalse\tGCA_958449795.1\t1"),
        "Main source of planted genome is not the donor:\n{}",
        genome_summary
    );
    assert!(output_dir.join("divergent_loci_upper_comparisons.tsv").exists());

    for locus in ["333279at4751", "96853at4751", "3487871at4751"] {
        let filename = format!("{}_masked.fasta.gz", locus);
        compare_fasta_files(
            &output_dir.join("masked_msas").join(&filename),
            &truth_dir.join(&filename),
        );
    }

    // plots are created for linear and log scales, with a panel per locus
    for name in ["divergent_loci_plots", "divergent_loci_plots_log", "divergent_loci_upper_plots", "divergent_loci_upper_plots_log"] {
        let png = fs::read(output_dir.join("plots").join(format!("{}.png", name))).unwrap();
        assert!(png.starts_with(b"\x89PNG"), "{}.png is not a PNG file", name);

        let svg = fs::read_to_string(output_dir.join("plots").join(format!("{}.svg", name))).unwrap();
        for locus in ["101409at4751", "333279at4751", "3487871at4751", "94226at4751", "96853at4751", "98953at4751"] {
            assert!(svg.contains(locus), "{}.svg has no panel for {}", name, locus);
        }
    }
}

// Comparisons with fewer aligned columns at the locus than --min-aligned-columns are excluded from both tests.
// The test set has 81 comparisons with 75-149 aligned columns, so a minimum of 150 removes them.
#[test]
fn test_mask_divergent_loci_min_aligned_columns() {
    let temp_dir = tempdir().unwrap();
    let output_dir = temp_dir.path();

    let mut cmd = cargo_bin_cmd!("gtrim");
    cmd.arg("mask-divergent-loci")
        .arg("-i")
        .arg("./tests/input_file/divergent_test")
        .arg("-o")
        .arg(output_dir)
        .arg("--window-size")
        .arg("40")
        .arg("--upper-window-size")
        .arg("40")
        .arg("--min-aligned-columns")
        .arg("150");

    let output = cmd.output().expect("failed to execute process");
    if !output.status.success() {
        println!("stderr: {}", String::from_utf8_lossy(&output.stderr));
        panic!("Command executed with failing error code");
    }

    for file in ["divergent_loci_comparisons.tsv", "divergent_loci_upper_comparisons.tsv"] {
        let content = fs::read_to_string(output_dir.join(file)).unwrap();
        let mut lines = content.lines();
        let col = lines.next().unwrap().split('\t').position(|c| c == "locus_aligned").unwrap();
        let aligned: Vec<u32> = lines.map(|l| l.split('\t').nth(col).unwrap().parse().unwrap()).collect();
        assert!(!aligned.is_empty(), "{} has no comparisons", file);
        assert!(aligned.iter().all(|&a| a >= 150), "{} has comparisons with fewer than 150 aligned columns", file);
    }
}
