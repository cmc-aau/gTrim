use std::fs::File;
use std::io::{BufRead, BufReader};

use anyhow::{Context, Result, bail};

/// A marker gene in a concatenated alignment.
#[derive(Debug, Clone, PartialEq)]
pub struct MarkerInfo {
    pub id: String,
    pub description: String,
    pub length: usize,
}

/// Read a marker info file (e.g. `gtdb_r<release>_bac120_markers_info.tsv`).
///
/// The file is tab-separated with a header containing the columns `Marker Id`, `Description`
/// and `Length (bp)` (matched case-insensitively, so GTDB-Tk's `Marker ID` also works). Markers
/// must be listed in the order they are concatenated in the alignment.
pub fn read_marker_info(path: &str) -> Result<Vec<MarkerInfo>> {
    let file = File::open(path).with_context(|| format!("Failed to open marker info file: {}", path))?;
    let mut lines = BufReader::new(file).lines();

    let header = match lines.next() {
        Some(line) => line.with_context(|| format!("Failed to read marker info file: {}", path))?,
        None => bail!("Marker info file is empty: {}", path),
    };
    let header: Vec<String> = header
        .trim_end_matches(['\r', '\n'])
        .split('\t')
        .map(|h| h.trim().to_lowercase())
        .collect();

    let column = |name: &str| {
        header
            .iter()
            .position(|h| h == name)
            .with_context(|| format!("Marker info file '{}' has no '{}' column", path, name))
    };
    let id_idx = column("marker id")?;
    let desc_idx = column("description")?;
    let len_idx = column("length (bp)")?;

    let mut markers = Vec::new();
    for (line_idx, line) in lines.enumerate() {
        let line = line.with_context(|| format!("Failed to read marker info file: {}", path))?;
        let line = line.trim_end_matches('\r');
        if line.trim().is_empty() {
            continue;
        }

        let line_num = line_idx + 2;
        let tokens: Vec<&str> = line.split('\t').collect();
        let field = |idx: usize| {
            tokens
                .get(idx)
                .map(|t| t.trim())
                .with_context(|| format!("Line {} of marker info file '{}' has too few columns", line_num, path))
        };

        let length_str = field(len_idx)?;
        let length = length_str.parse::<usize>().with_context(|| {
            format!("Invalid marker length '{}' on line {} of marker info file '{}'", length_str, line_num, path)
        })?;

        markers.push(MarkerInfo {
            id: field(id_idx)?.to_string(),
            description: field(desc_idx)?.to_string(),
            length,
        });
    }

    if markers.is_empty() {
        bail!("Marker info file contains no markers: {}", path);
    }

    Ok(markers)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn write_temp(content: &str) -> tempfile::NamedTempFile {
        let mut file = tempfile::NamedTempFile::new().unwrap();
        file.write_all(content.as_bytes()).unwrap();
        file
    }

    #[test]
    fn test_read_marker_info() {
        let file = write_temp(
            "Marker ID\tName\tDescription\tLength (bp)\n\
             PF00380.20\tRibosomal_S9\tRibosomal protein S9/S16\t121\n\
             TIGR00006\tTIGR00006\tS-adenosyl-methyltransferase MraW\t305\n\n",
        );
        let markers = read_marker_info(file.path().to_str().unwrap()).unwrap();
        assert_eq!(markers.len(), 2);
        assert_eq!(markers[0].id, "PF00380.20");
        assert_eq!(markers[0].description, "Ribosomal protein S9/S16");
        assert_eq!(markers[0].length, 121);
        assert_eq!(markers[1].length, 305);
    }

    #[test]
    fn test_read_marker_info_missing_column() {
        let file = write_temp("Marker Id\tName\tDescription\n");
        assert!(read_marker_info(file.path().to_str().unwrap()).is_err());
    }

    #[test]
    fn test_read_marker_info_bad_length() {
        let file = write_temp("Marker Id\tDescription\tLength (bp)\nM1\tdesc\tabc\n");
        assert!(read_marker_info(file.path().to_str().unwrap()).is_err());
    }
}
