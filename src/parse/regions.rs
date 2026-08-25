use mongodb::bson::Document;
use std::collections::HashMap;

const CODING_REGION_FLANK: i32 = 5_000;

#[derive(Debug, Clone)]
pub struct GenomicRegion {
    pub start: i32,
    pub end: i32,
    pub id: usize,
}

/// Builds merged genomic intervals around coding regions.
///
/// Each gene is expanded by 5 kb upstream and downstream. Overlapping
/// intervals on the same chromosome are merged into a single genomic region.
///
/// The genes are expected to have already been loaded for the requested
/// genome build. Only the values of `hgncid_to_gene` are used.
///
/// # Arguments
///
/// * `hgncid_to_gene` - Mapping from HGNC ID to gene documents.
///
/// # Returns
///
/// A map of chromosome names to merged genomic regions.
pub fn get_coding_intervals(
    hgncid_to_gene: &HashMap<i32, Document>,
) -> HashMap<String, Vec<GenomicRegion>> {
    let mut intervals: HashMap<String, Vec<GenomicRegion>> = HashMap::new();
    let mut next_region_id = 0;

    for gene in hgncid_to_gene.values() {
        let chromosome = match gene.get_str("chromosome") {
            Ok(chromosome) => chromosome.to_string(),
            Err(_) => continue,
        };

        let gene_start = match gene.get_i32("start") {
            Ok(start) => start,
            Err(_) => continue,
        };

        let gene_end = match gene.get_i32("end") {
            Ok(end) => end,
            Err(_) => continue,
        };

        let start = (gene_start - CODING_REGION_FLANK).max(1);
        let end = gene_end + CODING_REGION_FLANK;

        let chromosome_intervals = intervals.entry(chromosome).or_default();

        let mut merged_start = start;
        let mut merged_end = end;

        let mut overlapping = Vec::new();

        for (index, interval) in chromosome_intervals.iter().enumerate() {
            if interval.start < merged_end && interval.end > merged_start {
                merged_start = merged_start.min(interval.start);
                merged_end = merged_end.max(interval.end);
                overlapping.push(index);
            }
        }

        for index in overlapping.into_iter().rev() {
            chromosome_intervals.remove(index);
        }

        chromosome_intervals.push(GenomicRegion {
            start: merged_start,
            end: merged_end,
            id: next_region_id,
        });

        next_region_id += 1;
    }

    intervals
}

/// Finds the coding region overlapping a genomic interval.
///
/// The coding intervals are grouped by chromosome. The function returns the
/// identifier of the first coding region that overlaps the given interval,
/// or `None` if the interval does not overlap a coding region.
///
/// # Arguments
///
/// * `coding_intervals` - Coding regions grouped by chromosome.
/// * `chromosome` - Chromosome of the genomic interval.
/// * `start` - Start position of the genomic interval.
/// * `end` - End position of the genomic interval (exclusive).
///
/// # Returns
///
/// The identifier of the overlapping coding region, or `None` if the interval
/// is not within a coding region.
pub fn find_coding_region(
    coding_intervals: &HashMap<String, Vec<GenomicRegion>>,
    chromosome: &str,
    start: i32,
    end: i32,
) -> Option<usize> {
    coding_intervals
        .get(chromosome)
        .and_then(|regions| {
            regions
                .iter()
                .find(|region| start < region.end && end > region.start)
        })
        .map(|region| region.id)
}
