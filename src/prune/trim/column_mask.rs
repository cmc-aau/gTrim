use anyhow::{Context, Result, bail};

/// Format of a file indicating which columns of the full concatenated alignment were removed
/// before trimming.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MaskFormat {
    /// WitChi `*_pruned.tsv`: one row per removed column, with its 0-based position in the full
    /// alignment in the `Original Index` column.
    WitChi,
    /// A single line of `0`/`1` characters, one per column, where `1` = retained.
    Binary,
}

impl std::fmt::Display for MaskFormat {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            MaskFormat::WitChi => write!(f, "WitChi"),
            MaskFormat::Binary => write!(f, "binary mask"),
        }
    }
}

/// Columns of the full concatenated alignment retained before trimming.
#[derive(Debug, Clone, PartialEq)]
pub struct ColumnMask {
    pub format: MaskFormat,
    /// One entry per column of the full alignment; `true` = retained.
    pub retained: Vec<bool>,
}

/// Read a WitChi or binary mask file describing a full alignment of `full_len` columns.
///
/// The format is detected from the content: a single line of `0`s and `1`s is a binary mask, and a
/// tab-separated file with an `Original Index` column (matched case-insensitively) is a WitChi file.
pub fn read_column_mask(path: &str, full_len: usize) -> Result<ColumnMask> {
    let content = std::fs::read_to_string(path).with_context(|| format!("Failed to read mask file: {}", path))?;
    parse_column_mask(&content, path, full_len)
}

fn parse_column_mask(content: &str, path: &str, full_len: usize) -> Result<ColumnMask> {
    // (line number, line) of non-blank lines
    let mut lines = content
        .lines()
        .enumerate()
        .map(|(idx, line)| (idx + 1, line.trim_end_matches('\r')))
        .filter(|(_, line)| !line.trim().is_empty());

    let Some((_, first)) = lines.next() else {
        bail!("Mask file is empty: {}", path);
    };

    let first_trimmed = first.trim();
    if first_trimmed.bytes().all(|c| c == b'0' || c == b'1') {
        if let Some((line_num, _)) = lines.next() {
            bail!(
                "Binary mask file '{}' must contain a single line of 0s and 1s, but line {} is not empty.",
                path,
                line_num
            );
        }

        let retained: Vec<bool> = first_trimmed.bytes().map(|c| c == b'1').collect();
        if retained.len() != full_len {
            bail!(
                "Binary mask file '{}' has {} columns, but the markers in the marker info file have a summed length of {} columns.",
                path,
                retained.len(),
                full_len
            );
        }

        return Ok(ColumnMask {
            format: MaskFormat::Binary,
            retained,
        });
    }

    let index_col = first
        .split('\t')
        .position(|h| h.trim().eq_ignore_ascii_case("original index"))
        .with_context(|| {
            format!(
                "Mask file '{}' is neither a binary mask (a single line of 0s and 1s) nor a WitChi file (a tab-separated file with an 'Original Index' column).",
                path
            )
        })?;

    let mut retained = vec![true; full_len];
    for (line_num, line) in lines {
        let token = line
            .split('\t')
            .nth(index_col)
            .map(str::trim)
            .with_context(|| format!("Line {} of WitChi file '{}' has too few columns", line_num, path))?;
        let col_idx = token.parse::<usize>().with_context(|| {
            format!("Invalid column index '{}' on line {} of WitChi file '{}'", token, line_num, path)
        })?;

        if col_idx >= full_len {
            bail!(
                "Column index {} on line {} of WitChi file '{}' is outside the full alignment: indices are 0-based and the markers in the marker info file have a summed length of {} columns.",
                col_idx,
                line_num,
                path,
                full_len
            );
        }
        if !retained[col_idx] {
            bail!("Column index {} is listed more than once in WitChi file '{}' (line {}).", col_idx, path, line_num);
        }
        retained[col_idx] = false;
    }

    Ok(ColumnMask {
        format: MaskFormat::WitChi,
        retained,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    const WITCHI_HEADER: &str = "Iteration\tOriginal Index\tGlobal ChiScore\tInitial Significant ChiScore\tChiScore Difference\tSignificant taxa\tNull-Pvalue\tAlgorithm\n";

    fn parse(content: &str, full_len: usize) -> Result<ColumnMask> {
        parse_column_mask(content, "test", full_len)
    }

    #[test]
    fn test_binary_mask() {
        let mask = parse("10110\n", 5).unwrap();
        assert_eq!(mask.format, MaskFormat::Binary);
        assert_eq!(mask.retained, vec![true, false, true, true, false]);

        // Windows line ending, no trailing newline
        assert_eq!(parse("01\r\n", 2).unwrap().retained, vec![false, true]);
        assert_eq!(parse("11", 2).unwrap().retained, vec![true, true]);
    }

    #[test]
    fn test_binary_mask_wrong_length() {
        let err = parse("1011\n", 5).unwrap_err().to_string();
        assert!(err.contains("has 4 columns"), "{}", err);
    }

    #[test]
    fn test_binary_mask_multiple_lines() {
        assert!(parse("101\n010\n", 3).is_err());
    }

    #[test]
    fn test_witchi() {
        let content = format!("{}1\t3\t10.0\t9.0\t1.0\t5\t\tquartic\n2\t0\t9.0\t8.0\t1.0\t4\t\tquartic\n", WITCHI_HEADER);
        let mask = parse(&content, 5).unwrap();
        assert_eq!(mask.format, MaskFormat::WitChi);
        assert_eq!(mask.retained, vec![false, true, true, false, true]);
    }

    #[test]
    fn test_witchi_no_rows() {
        let mask = parse(WITCHI_HEADER, 3).unwrap();
        assert_eq!(mask.retained, vec![true; 3]);
    }

    #[test]
    fn test_witchi_index_out_of_range() {
        let content = format!("{}1\t5\t10.0\t9.0\t1.0\t5\t\tquartic\n", WITCHI_HEADER);
        let err = parse(&content, 5).unwrap_err().to_string();
        assert!(err.contains("outside the full alignment"), "{}", err);
    }

    #[test]
    fn test_witchi_duplicate_index() {
        let content = format!("{}1\t2\tx\n2\t2\tx\n", WITCHI_HEADER);
        let err = parse(&content, 5).unwrap_err().to_string();
        assert!(err.contains("more than once"), "{}", err);
    }

    #[test]
    fn test_witchi_bad_index() {
        let content = format!("{}1\tabc\tx\n", WITCHI_HEADER);
        assert!(parse(&content, 5).is_err());
    }

    #[test]
    fn test_unknown_format() {
        let err = parse("Column\tStatus\n1\tremoved\n", 5).unwrap_err().to_string();
        assert!(err.contains("neither a binary mask"), "{}", err);
        assert!(parse("", 5).is_err());
    }
}
