use std::collections::{BTreeSet, HashMap};
use std::path::{Path, PathBuf};

use anyhow::{Context, Result};
use log::info;
use ndarray::Array2;

use crate::alignment_reader::AlignmentReader;
use crate::alignment_writer::strip_seq_extensions;
use crate::sequence_type_detector::SequenceTypeDetector;

/// Encoded value for gaps, ambiguous residues, and missing loci.
pub const GAP_CODE: u8 = 0;

/// A single locus (gene) MSA within the concatenated genome matrix.
pub struct Locus {
    pub name: String,
    pub file: PathBuf,
    /// First column of this locus in the concatenated matrix.
    pub offset: usize,
    /// Number of alignment columns.
    pub len: usize,
}

/// Collection of per-locus MSAs, concatenated into a genome x column matrix.
pub struct LocusSet {
    pub loci: Vec<Locus>,
    pub genome_ids: Vec<String>,
    /// Genomes x concatenated columns. Valid residues are stored as their uppercase
    /// byte and gaps, ambiguous residues, and missing loci as GAP_CODE.
    pub matrix: Array2<u8>,
    /// Genomes x loci count of valid (non-gap) residues. Zero if the locus is absent.
    pub nongap: Array2<u32>,
    /// Genomes x loci indicating if the genome has a sequence for the locus.
    pub present: Array2<bool>,
    /// Total number of valid residues across all loci for each genome.
    pub nongap_total: Vec<u32>,
    pub is_dna: bool,
}

/// Get genome ID from a FASTA header (first whitespace delimited token).
pub fn genome_id_from_header(header: &str) -> &str {
    header.split_whitespace().next().unwrap_or("")
}

/// Get sorted list of MSA files in a directory.
pub fn msa_files(msa_dir: &str) -> Result<Vec<PathBuf>> {
    let mut files = Vec::new();
    for entry in std::fs::read_dir(msa_dir)
        .with_context(|| format!("Failed to read MSA directory: {}", msa_dir))?
    {
        let path = entry?.path();
        if !path.is_file() {
            continue;
        }

        let filename = path.file_name().and_then(|n| n.to_str()).unwrap_or("");
        let stem = filename.trim_end_matches(".gz");
        if stem.ends_with(".faa") || stem.ends_with(".fna") || stem.ends_with(".fasta") {
            files.push(path);
        }
    }

    files.sort();
    Ok(files)
}

impl LocusSet {
    /// Read all MSAs in a directory.
    pub fn from_dir(msa_dir: &str) -> Result<Self> {
        let files = msa_files(msa_dir)?;
        if files.is_empty() {
            anyhow::bail!("No MSA files (.faa, .fna, .fasta, optionally .gz) found in: {}", msa_dir);
        }

        // read and encode each locus
        let mut loci = Vec::with_capacity(files.len());
        let mut encoded_loci = Vec::with_capacity(files.len());
        let mut all_genome_ids = BTreeSet::new();
        let mut is_dna = None;
        let mut offset = 0;
        for file in files {
            let (locus_genome_ids, encoded, locus_is_dna) = Self::read_locus(&file)?;

            if *is_dna.get_or_insert(locus_is_dna) != locus_is_dna {
                anyhow::bail!("MSA {} has a different sequence type to previous MSAs.", file.display());
            }

            all_genome_ids.extend(locus_genome_ids.iter().cloned());

            let filename = file.file_name().and_then(|n| n.to_str()).unwrap_or("");
            let len = encoded.ncols();
            loci.push(Locus {
                name: strip_seq_extensions(filename).to_string(),
                file,
                offset,
                len,
            });
            encoded_loci.push((locus_genome_ids, encoded));
            offset += len;
        }

        // build concatenated genome x column matrix
        let genome_ids: Vec<String> = all_genome_ids.into_iter().collect();
        let genome_idx: HashMap<&str, usize> = genome_ids
            .iter()
            .enumerate()
            .map(|(idx, id)| (id.as_str(), idx))
            .collect();

        let mut matrix = Array2::from_elem((genome_ids.len(), offset), GAP_CODE);
        let mut nongap = Array2::zeros((genome_ids.len(), loci.len()));
        let mut present = Array2::from_elem((genome_ids.len(), loci.len()), false);
        for (locus_idx, (locus_genome_ids, encoded)) in encoded_loci.into_iter().enumerate() {
            let locus = &loci[locus_idx];
            for (row_idx, genome_id) in locus_genome_ids.iter().enumerate() {
                let g = genome_idx[genome_id.as_str()];
                let row = encoded.row(row_idx);
                for (col_idx, &code) in row.iter().enumerate() {
                    matrix[[g, locus.offset + col_idx]] = code;
                }
                nongap[[g, locus_idx]] = row.iter().filter(|&&c| c != GAP_CODE).count() as u32;
                present[[g, locus_idx]] = true;
            }
        }

        let nongap_total = nongap.rows().into_iter().map(|r| r.sum()).collect();

        info!(
            "Read {} loci spanning {} genomes and {} concatenated columns.",
            loci.len(),
            genome_ids.len(),
            offset
        );

        Ok(Self {
            loci,
            genome_ids,
            matrix,
            nongap,
            present,
            nongap_total,
            is_dna: is_dna.unwrap_or(false),
        })
    }

    /// Read a single locus MSA and encode residues. Returns genome IDs, encoded MSA, and if it is DNA.
    fn read_locus(file: &Path) -> Result<(Vec<String>, Array2<u8>, bool)> {
        let file_str = file.to_string_lossy().to_string();
        let reader = AlignmentReader::new(file_str.clone());
        let (seq_ids, alignment_array) = reader.run()?;

        let (is_dna, char_set) = SequenceTypeDetector::detect(&alignment_array)
            .with_context(|| format!("Failed to detect sequence type of: {}", file_str))?;

        let mut valid = [false; 256];
        for ch in char_set {
            valid[ch as usize] = true;
        }

        let encoded = alignment_array.mapv(|ch| {
            if ch.is_ascii() && valid[ch as usize] {
                ch as u8
            } else {
                GAP_CODE
            }
        });

        let mut genome_ids = Vec::with_capacity(seq_ids.len());
        let mut seen = BTreeSet::new();
        for header in &seq_ids {
            let genome_id = genome_id_from_header(header).to_string();
            if !seen.insert(genome_id.clone()) {
                anyhow::bail!("Genome {} occurs multiple times in: {}", genome_id, file_str);
            }
            genome_ids.push(genome_id);
        }

        Ok((genome_ids, encoded, is_dna))
    }

    pub fn num_genomes(&self) -> usize {
        self.genome_ids.len()
    }

    /// Encoded residues of a genome across all loci.
    pub fn genome_slice(&self, genome_idx: usize) -> &[u8] {
        self.matrix.row(genome_idx).to_slice().expect("matrix is in standard layout")
    }

    /// Encoded residues of a genome for a single locus.
    pub fn locus_slice(&self, genome_idx: usize, locus_idx: usize) -> &[u8] {
        let locus = &self.loci[locus_idx];
        &self.genome_slice(genome_idx)[locus.offset..locus.offset + locus.len]
    }
}
