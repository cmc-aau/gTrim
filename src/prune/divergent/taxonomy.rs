use std::collections::HashMap;

use anyhow::{Context, Result};

/// Taxonomic classification of genomes, read from a tab-separated file of genome ID and taxonomy string. The
/// taxonomy string lists ranks from highest to lowest separated by semicolons, e.g. GTDB format:
/// `d__Eukaryota;p__Ascomycota;c__Sordariomycetes;o__Hypocreales;f__Nectriaceae;g__Fusarium;s__Fusarium oxysporum`
pub struct Taxonomy {
    taxa: HashMap<String, Vec<String>>,
}

/// A rank without a name (e.g. `g__`) is unclassified at that rank.
fn is_named(rank: &str) -> bool {
    !rank.is_empty() && !rank.ends_with("__")
}

/// Remove the GTDB database prefix (RS_ or GB_) from a genome ID.
fn strip_db_prefix(genome_id: &str) -> &str {
    genome_id
        .strip_prefix("RS_")
        .or_else(|| genome_id.strip_prefix("GB_"))
        .unwrap_or(genome_id)
}

impl Taxonomy {
    /// Read taxonomy file. Empty lines and lines starting with '#' are ignored.
    pub fn from_file(path: &str) -> Result<Self> {
        let content = std::fs::read_to_string(path).with_context(|| format!("Failed to read taxonomy file: {}", path))?;
        Self::parse(&content).with_context(|| format!("Failed to parse taxonomy file: {}", path))
    }

    fn parse(content: &str) -> Result<Self> {
        let mut taxa = HashMap::new();
        for (line_idx, line) in content.lines().enumerate() {
            let line = line.trim_end();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }

            let (genome_id, taxonomy) = line
                .split_once('\t')
                .with_context(|| format!("Line {} is not <genome ID><tab><taxonomy>", line_idx + 1))?;
            let ranks = taxonomy.split(';').map(|r| r.trim().to_string()).collect();
            taxa.insert(strip_db_prefix(genome_id.trim()).to_string(), ranks);
        }

        Ok(Self { taxa })
    }

    /// Check if a genome has a taxonomy.
    pub fn contains(&self, genome_id: &str) -> bool {
        self.taxa.contains_key(genome_id)
    }

    /// Lowest taxon shared by all genomes. The shared lineage ends at the first rank without a name. None if
    /// there are no genomes, any genome has no taxonomy, or the genomes share no named taxon.
    pub fn lowest_common_taxon<'a>(&self, genome_ids: impl IntoIterator<Item = &'a str>) -> Option<&str> {
        let mut lineages = genome_ids.into_iter().map(|id| self.taxa.get(id));
        let first = lineages.next()??;

        let mut depth = first.len();
        for lineage in lineages {
            let lineage = lineage?;
            depth = depth.min(first.iter().zip(lineage).take_while(|(a, b)| a == b).count());
        }
        let depth = first[..depth].iter().position(|r| !is_named(r)).unwrap_or(depth);

        depth.checked_sub(1).map(|d| first[d].as_str())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn taxonomy() -> Taxonomy {
        Taxonomy::parse(
            "# genome\ttaxonomy\n\
             GB_GCA_1.1\td__E;p__A;c__S;o__H;f__N;g__F;s__F x\n\
             RS_GCF_2.1\td__E;p__A;c__S;o__H;f__N;g__F;s__F y\n\
             GCA_3.1\td__E;p__A;c__S;o__G;f__;g__;s__\n\
             GCA_4.1\td__E;p__B;c__T;o__U;f__V;g__W;s__W z\n",
        )
        .unwrap()
    }

    #[test]
    fn test_lowest_common_taxon() {
        let t = taxonomy();
        assert_eq!(t.lowest_common_taxon(["GCA_1.1"]), Some("s__F x"));
        assert_eq!(t.lowest_common_taxon(["GCA_1.1", "GCF_2.1"]), Some("g__F"));
        assert_eq!(t.lowest_common_taxon(["GCA_1.1", "GCF_2.1", "GCA_3.1"]), Some("c__S"));
        assert_eq!(t.lowest_common_taxon(["GCA_1.1", "GCA_4.1"]), Some("d__E"));
        assert_eq!(t.lowest_common_taxon(["GCA_1.1", "GCA_missing"]), None);
        assert_eq!(t.lowest_common_taxon([]), None);
    }

    #[test]
    fn test_unnamed_ranks_end_lineage() {
        // GCA_3.1 is unclassified below order, so f__ isn't a shared taxon
        let t = taxonomy();
        assert_eq!(t.lowest_common_taxon(["GCA_3.1"]), Some("o__G"));
        assert!(Taxonomy::parse("GCA_1.1 no tab\n").is_err());
    }
}
