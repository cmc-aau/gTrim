use std::fs::File;
use std::io::{Write, BufWriter};

use anyhow::{Context, Result};
use ndarray::Array2;
use flate2::write::GzEncoder;
use flate2::Compression;


/// Check if input file is compressed.
pub fn check_compressed_input(file: &str) -> bool {
    file.ends_with(".gz")
}

/// Get output prefix which consists of the output path and the input filename without any extensions.
pub fn get_output_prefix(input_file: &str, out_dir: &str) -> Result<String> {
    // get the filename from the input path
    let mut filename = std::path::Path::new(input_file)
        .file_name()
        .and_then(|name| name.to_str())
        .ok_or_else(|| anyhow::anyhow!("Input path '{}' does not contain a valid filename", input_file))?;

    // create owned stem and strip known extensions
    filename = filename
        .trim_end_matches(".gz")
        .trim_end_matches(".faa")
        .trim_end_matches(".fna")
        .trim_end_matches(".fasta");

    // join output directory with the sanitized stem and return as String
    let base_path = std::path::Path::new(out_dir)
        .join(filename)
        .to_string_lossy()
        .into_owned();

    Ok(base_path)
}

/// Write filtered alignment to FASTA file.
pub fn write_filtered_alignment(
    seq_ids: &[String],
    alignment_array: &Array2<char>,
    mask: &[bool],
    out_file: &str,
    compress: bool,
) -> Result<()> {
    // Append .gz extension if compression is requested
    let out_path = if compress {
        format!("{}.gz", out_file)
    } else {
        out_file.to_string()
    };

    let file = File::create(&out_path)
        .with_context(|| format!("Failed to create file: {}", out_path))?;

    let mut writer: Box<dyn Write> = if compress {
        Box::new(GzEncoder::new(file, Compression::default()))
    } else {
        Box::new(BufWriter::new(file))
    };

    for (row_idx, id) in seq_ids.iter().enumerate() {
        writeln!(writer, ">{}", id)?;
        let seq: String = alignment_array
            .row(row_idx)
            .iter()
            .copied()
            .enumerate()
            .filter_map(|(i, ch)| if mask[i] { Some(ch) } else { None })
            .collect();
        writeln!(writer, "{}", seq)?;
    }

    writer.flush()?;

    Ok(())
}

/// Write out mask (true = 1; false = 0).
pub fn write_mask(mask: &[bool], out_file: &str) -> Result<()> {
    // Convert bool mask to a string of '1' and '0' characters
    let mask_string: String = mask.iter().map(|&b| if b { '1' } else { '0' }).collect();

    let mut file = File::create(out_file).with_context(|| format!("Failed to create file: {}", out_file))?;
    writeln!(file, "{}", mask_string)?;

    Ok(())
}