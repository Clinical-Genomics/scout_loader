use mongodb::bson::{Document, doc};
use scout_loader::parse::regions::{GenomicRegion, find_coding_region, get_coding_intervals};
use std::collections::HashMap;

fn gene(chromosome: &str, start: i32, end: i32) -> Document {
    doc! { "chromosome": chromosome, "start": start, "end": end }
}

fn genes(genes: Vec<Document>) -> HashMap<i32, Document> {
    genes
        .into_iter()
        .enumerate()
        .map(|(hgnc_id, gene)| (hgnc_id as i32, gene))
        .collect()
}

/// Region bounds per chromosome, sorted, ignoring the (order-dependent) ids.
fn bounds(intervals: &HashMap<String, Vec<GenomicRegion>>, chromosome: &str) -> Vec<(i32, i32)> {
    let mut bounds: Vec<(i32, i32)> = intervals[chromosome]
        .iter()
        .map(|region| (region.start, region.end))
        .collect();
    bounds.sort();
    bounds
}

#[test]
fn test_get_coding_intervals_adds_flanks() {
    // GIVEN a single gene
    let intervals = get_coding_intervals(&genes(vec![gene("1", 100_000, 200_000)]));

    // THEN the region is padded by 5 kb on each side
    assert_eq!(bounds(&intervals, "1"), [(95_000, 205_000)]);
}

#[test]
fn test_get_coding_intervals_clamps_start() {
    // GIVEN a gene close to the start of the chromosome
    let intervals = get_coding_intervals(&genes(vec![gene("MT", 577, 647)]));

    // THEN the region starts at position 1
    assert_eq!(bounds(&intervals, "MT"), [(1, 5_647)]);
}

#[test]
fn test_get_coding_intervals_merges_overlapping_genes() {
    // GIVEN two genes whose flanked regions overlap, one far away,
    // and one on another chromosome
    let intervals = get_coding_intervals(&genes(vec![
        gene("1", 100_000, 110_000),
        gene("1", 118_000, 130_000),
        gene("1", 500_000, 510_000),
        gene("2", 100_000, 110_000),
    ]));

    // THEN the overlapping regions are merged, per chromosome
    assert_eq!(
        bounds(&intervals, "1"),
        [(95_000, 135_000), (495_000, 515_000)]
    );
    assert_eq!(bounds(&intervals, "2"), [(95_000, 115_000)]);
}

#[test]
fn test_get_coding_intervals_merges_bridging_gene() {
    // GIVEN two separate genes and a third gene bridging them
    let intervals = get_coding_intervals(&genes(vec![
        gene("1", 100_000, 110_000),
        gene("1", 130_000, 140_000),
        gene("1", 110_000, 130_000),
    ]));

    // THEN everything collapses into one region, whatever the iteration order
    assert_eq!(bounds(&intervals, "1"), [(95_000, 145_000)]);
}

#[test]
fn test_get_coding_intervals_skips_incomplete_genes() {
    // GIVEN genes missing chromosome, start or end (or with wrong types)
    let intervals = get_coding_intervals(&genes(vec![
        doc! { "start": 1_000, "end": 2_000 },
        doc! { "chromosome": "1", "end": 2_000 },
        doc! { "chromosome": "1", "start": 1_000 },
        doc! { "chromosome": "1", "start": "1000", "end": 2_000 },
    ]));

    // THEN no intervals are created
    assert!(intervals.is_empty());
}

#[test]
fn test_find_coding_region() {
    // GIVEN one region on chromosome 1
    let intervals = HashMap::from([(
        "1".to_string(),
        vec![GenomicRegion {
            start: 1_000,
            end: 2_000,
            id: 7,
        }],
    )]);

    // THEN overlapping intervals return the region id
    assert_eq!(find_coding_region(&intervals, "1", 1_500, 1_501), Some(7));
    assert_eq!(find_coding_region(&intervals, "1", 500, 1_001), Some(7));
    assert_eq!(find_coding_region(&intervals, "1", 1_999, 3_000), Some(7));

    // AND end is exclusive on both sides
    assert_eq!(find_coding_region(&intervals, "1", 500, 1_000), None);
    assert_eq!(find_coding_region(&intervals, "1", 2_000, 2_001), None);

    // AND other chromosomes do not match
    assert_eq!(find_coding_region(&intervals, "2", 1_500, 1_501), None);
}
