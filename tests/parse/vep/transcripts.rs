//! VEP transcript parsing tests, ported from scout's
//! `tests/parse/test_parse_transcripts.py`.
//!
//! The CSQ fixtures are copied from scout's `tests/parse/conftest.py`
//! (`vep_csq_header`/`vep_csq`) and `test_parse_transcripts.py` (ABCC2).

use crate::vcf_record::vep_record;
use mongodb::bson::{Bson, Document, doc};
use scout_loader::parse::vep::transcripts::parse_vep_transcripts;

/// scout `vep_csq_header`: SCN5A, VEP with MANE, 1000G, gnomAD and gnomAD-mt fields.
const SCN5A_CSQ_FORMAT: &str = "Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|VARIANT_CLASS|SYMBOL_SOURCE|HGNC_ID|CANONICAL|MANE_SELECT|MANE_PLUS_CLINICAL|TSL|APPRIS|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|UNIPROT_ISOFORM|REFSEQ_MATCH|SOURCE|REFSEQ_OFFSET|GIVEN_REF|USED_REF|BAM_EDIT|GENE_PHENO|SIFT|PolyPhen|DOMAINS|miRNA|HGVS_OFFSET|AF|AFR_AF|AMR_AF|EAS_AF|EUR_AF|SAS_AF|AA_AF|EA_AF|gnomAD_AF|gnomAD_AFR_AF|gnomAD_AMR_AF|gnomAD_ASJ_AF|gnomAD_EAS_AF|gnomAD_FIN_AF|gnomAD_NFE_AF|gnomAD_OTH_AF|gnomAD_SAS_AF|MAX_AF|MAX_AF_POPS|CLIN_SIG|SOMATIC|PHENO|PUBMED|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|TRANSCRIPTION_FACTORS|CADD_PHRED|CADD_RAW|LoFtool|MES-NCSS_downstream_acceptor|MES-NCSS_downstream_donor|MES-NCSS_upstream_acceptor|MES-NCSS_upstream_donor|MES-SWA_acceptor_alt|MES-SWA_acceptor_diff|MES-SWA_acceptor_ref|MES-SWA_acceptor_ref_comp|MES-SWA_donor_alt|MES-SWA_donor_diff|MES-SWA_donor_ref|MES-SWA_donor_ref_comp|MaxEntScan_alt|MaxEntScan_diff|MaxEntScan_ref|gnomAD_mt_AF_hom|gnomAD_mt_AF_het|phyloP100way|phastCons";

/// scout `vep_csq`. Has both MANE_SELECT and MANE_PLUS_CLINICAL set for testing.
const SCN5A_CSQ: &str = "C|synonymous_variant|LOW|SCN5A|6331|Transcript|NM_000335.5|protein_coding|17/28||NM_000335.5:c.3183A>G|NP_000326.2:p.Glu1061%3D|3392|3183|1061|E|gaA/gaG|rs7430407&CS1510231&COSV61125864||-1||SNV|EntrezGene|HGNC:10593|YES|NM_000335.5|NM_001099404.2||||NP_000326.2||||||RefSeq||T|T||||||||0.9242|0.8661|0.9164|0.997|0.9085|0.9438|0.8889|0.883|0.9139|0.8876|0.9444|0.8688|0.9993|0.9275|0.8895|0.9041|0.9398|0.9993|gnomAD_EAS|likely_benign&benign|0&0&1|1&1&1|25741868&24033266&27182706&26401487&25401102&20004937||||||1.179|-0.072561|0.000413|5.050|8.066|9.991|10.355|-11.882|-0.089|-11.970|-11.970|-1.050|2.025|0.976|0.976||||0.2|0.0123|-1.08899998664856|0.0209999997168779";

/// ABCC2 missense variant, three transcripts (Ensembl + two RefSeq).
const ABCC2_CSQ_FORMAT: &str = "Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|SYMBOL_SOURCE|HGNC_ID|CANONICAL|TSL|APPRIS|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|REFSEQ_MATCH|SOURCE|GIVEN_REF|USED_REF|BAM_EDIT|SIFT|PolyPhen|DOMAINS|HGVS_OFFSET|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|MES-NCSS_downstream_acceptor|MES-NCSS_downstream_donor|MES-NCSS_upstream_acceptor|MES-NCSS_upstream_donor|MES-SWA_acceptor_alt|MES-SWA_acceptor_diff|MES-SWA_acceptor_ref|MES-SWA_acceptor_ref_comp|MES-SWA_donor_alt|MES-SWA_donor_diff|MES-SWA_donor_ref|MES-SWA_donor_ref_comp|MaxEntScan_alt|MaxEntScan_diff|MaxEntScan_ref|GERP++_NR|GERP++_RS|REVEL_rankscore|phastCons100way_vertebrate|phyloP100way_vertebrate|LoFtool|ExACpLI|CLINVAR|CLINVAR_CLNSIG|CLINVAR_CLNVID|CLINVAR_CLNREVSTAT|genomic_superdups_frac_match";

const ABCC2_CSQ: &str = "A|missense_variant|MODERATE|ABCC2|ENSG00000023839|Transcript|ENST00000370449|protein_coding|10/32||ENST00000370449.4:c.1249G>A|ENSP00000359478.4:p.Val417Ile|1362|1249|417|V/I|Gtt/Att|||1||HGNC|53|YES|||CCDS7484.1|ENSP00000359478|Q92887||UPI000013D6CA||Ensembl|G|G||tolerated|benign|PROSITE_profiles:PS50929&hmmpanther:PTHR24223&hmmpanther:PTHR24223:SF176&TIGRFAM_domain:TIGR00957&Pfam_domain:PF00664&Gene3D:2hydA01&Superfamily_domains:SSF90123||||||12.834|9.654|8.373|10.218|-7.208|1.111|-6.097|-6.097|-12.417|-1.215|-9.514|-13.633||||5.48|-10.3|0.54096|0.000000|-0.327000|0.0978|0.00|255992|Benign|255992|criteria_provided&_multiple_submitters&_no_conflicts|,A|missense_variant|MODERATE|ABCC2|1244|Transcript|NM_000392.3|protein_coding|10/32||NM_000392.3:c.1249G>A|NP_000383.1:p.Val417Ile|1388|1249|417|V/I|Gtt/Att|||1||EntrezGene|53|YES||||NP_000383.1|||||RefSeq|G|G||tolerated|benign|||||||12.834|9.654|8.373|10.218|-7.208|1.111|-6.097|-6.097|-12.417|-1.215|-9.514|-13.633||||5.48|-10.3|0.54096|0.000000|-0.327000|0.0978|0.00|255992|Benign|255992|criteria_provided&_multiple_submitters&_no_conflicts|,A|missense_variant|MODERATE|ABCC2|1244|Transcript|NM_000392.4|protein_coding|10/32||NM_000392.4:c.1249G>A|NP_000383.1:p.Val417Ile|1496|1249|417|V/I|Gtt/Att|||1||EntrezGene|53|YES||||NP_000383.1||||rseq_mrna_nonmatch&rseq_cds_mismatch|RefSeq|G|G|OK|tolerated|benign|||||||12.834|9.654|8.373|10.218|-7.208|1.111|-6.097|-6.097|-12.417|-1.215|-9.514|-13.633||||5.48|-10.3|0.54096|0.000000|-0.327000|0.0978|0.00|255992|Benign|255992|criteria_provided&_multiple_submitters&_no_conflicts|";

/// Parse the transcripts of a single-record VCF with the given CSQ.
fn transcripts(csq_format: &str, csq: &str) -> Vec<Document> {
    transcripts_with_variant(csq_format, csq, &[], "").0
}

fn transcripts_with_variant(
    csq_format: &str,
    csq: &str,
    extra_meta: &[&str],
    extra_info: &str,
) -> (Vec<Document>, Document) {
    let (record, field_indices) = vep_record(csq_format, csq, extra_meta, extra_info);
    let mut variant = Document::new();
    let transcripts = parse_vep_transcripts(&record, &field_indices, &mut variant);

    (transcripts, variant)
}

fn strings(values: &[&str]) -> Bson {
    Bson::Array(values.iter().map(|value| Bson::from(*value)).collect())
}

#[test]
fn test_parse_cosmic_csq_cosmic() {
    // GIVEN a CSQ entry containing a dedicated COSMIC field
    let transcripts = transcripts(
        "Feature|Consequence|COSMIC|COSMIC_CDS|COSMIC_GENE|COSMIC_STRAND|COSMIC_CNT|COSMIC_AA",
        "ENST00000264233|missense_variant|COSV99072206|c.913T>G|POLQ||2|p.S305A",
    );

    // THEN the COSMIC id is parsed
    assert_eq!(transcripts.len(), 1);
    assert_eq!(
        transcripts[0].get("cosmic"),
        Some(&strings(&["COSV99072206"]))
    );
}

#[test]
fn test_parse_transcripts() {
    // GIVEN three ABCC2 transcripts
    let transcripts = transcripts(ABCC2_CSQ_FORMAT, ABCC2_CSQ);

    // THEN all are parsed, in order
    let ids: Vec<&str> = transcripts
        .iter()
        .map(|transcript| transcript.get_str("transcript_id").unwrap())
        .collect();
    assert_eq!(ids, ["ENST00000370449", "NM_000392.3", "NM_000392.4"]);

    // AND the Ensembl transcript has the expected annotations
    let transcript = &transcripts[0];
    assert_eq!(transcript.get_i32("hgnc_id"), Ok(53));
    assert_eq!(transcript.get_str("protein_id"), Ok("ENSP00000359478"));
    assert_eq!(transcript.get_str("sift_prediction"), Ok("tolerated"));
    assert_eq!(transcript.get_str("polyphen_prediction"), Ok("benign"));
    assert_eq!(
        transcript.get("functional_annotations"),
        Some(&strings(&["missense_variant"]))
    );
    assert_eq!(
        transcript.get("region_annotations"),
        Some(&strings(&["exonic"]))
    );
    assert_eq!(transcript.get_str("gerp"), Ok("-10.3"));
    assert_eq!(transcript.get_str("phast"), Ok("0.000000"));
    assert_eq!(transcript.get_str("phylop"), Ok("-0.327000"));
    assert_eq!(transcript.get_f64("revel_rankscore"), Ok(0.54096));
    assert_eq!(transcript.get_str("swiss_prot"), Ok("Q92887"));
    assert_eq!(transcript.get_str("coding_sequence_name"), Ok("c.1249G>A"));
    assert_eq!(
        transcript.get_str("protein_sequence_name"),
        Ok("p.Val417Ile")
    );
    assert_eq!(transcript.get_str("biotype"), Ok("protein_coding"));
    assert_eq!(transcript.get_str("exon"), Ok("10/32"));
    assert_eq!(transcript.get_str("strand"), Ok("+"));
    assert_eq!(transcript.get_bool("is_canonical"), Ok(true));
    assert_eq!(transcript.get_str("clinvar_clnvid"), Ok("255992"));
    assert_eq!(transcript.get_str("clinvar_clnsig"), Ok("benign"));
    assert_eq!(
        transcript.get_str("clinvar_revstat"),
        Ok("criteria_provided&_multiple_submitters&_no_conflicts")
    );

    // AND the domains are parsed; the last hmmpanther entry wins
    assert_eq!(transcript.get_str("prosite_profile"), Ok("PS50929"));
    assert_eq!(transcript.get_str("pfam_domain"), Ok("PF00664"));
    assert_eq!(transcript.get_str("panther_domain"), Ok("PTHR24223:SF176"));

    // AND a transcript without SWISSPROT gets "unknown"
    assert_eq!(transcripts[1].get_str("swiss_prot"), Ok("unknown"));
}

#[test]
fn test_parse_functional_annotation_and_hgnc() {
    // GIVEN the scout vep_csq fixture (HGNC id with "HGNC:" prefix)
    let transcripts = transcripts(SCN5A_CSQ_FORMAT, SCN5A_CSQ);
    assert_eq!(transcripts.len(), 1);

    // THEN functional annotation and HGNC id are parsed
    for transcript in &transcripts {
        assert_eq!(
            transcript.get("functional_annotations"),
            Some(&strings(&["synonymous_variant"]))
        );
        assert_eq!(transcript.get_i32("hgnc_id"), Ok(10593));
        assert_eq!(transcript.get_str("hgnc_symbol"), Ok("SCN5A"));
        assert_eq!(transcript.get_str("strand"), Ok("-"));
        assert_eq!(
            transcript.get_str("protein_sequence_name"),
            Ok("p.Glu1061%3D")
        );
    }
}

#[test]
fn test_parse_vep_freq_thousand_g() {
    // GIVEN a transcript with the 1000GAF frequency
    let transcripts = transcripts(
        "Allele|Consequence|Feature|1000GAF",
        "C|missense_variant|ENST1|0.01",
    );

    // THEN it is parsed as thousand_g_maf
    assert_eq!(transcripts[0].get_f64("thousand_g_maf"), Ok(0.01));
}

#[test]
fn test_parse_vep_freq_scn5a() {
    // GIVEN the scout vep_csq fixture with AF, gnomAD_AF, population AFs and gnomAD-mt
    let transcripts = transcripts(SCN5A_CSQ_FORMAT, SCN5A_CSQ);
    let transcript = &transcripts[0];

    // THEN the 1000G AF is parsed (test_parse_vep_freq_thousand_g_alt)
    assert_eq!(transcript.get_f64("thousand_g_maf"), Ok(0.9242));

    // AND the gnomAD AF and max are parsed (test_parse_vep_freq_gnomad)
    assert_eq!(transcript.get_f64("gnomad_maf"), Ok(0.9139));
    assert_eq!(transcript.get_f64("gnomad_max"), Ok(0.9993));

    // AND the max over the remaining population frequencies (incl. MAX_AF)
    assert_eq!(transcript.get_f64("thousandg_max"), Ok(0.9993));

    // AND the mitochondrial gnomAD frequencies (test_parse_vep_freq_mtgnomad)
    assert_eq!(transcript.get_f64("gnomad_mt_homoplasmic"), Ok(0.2));
    assert_eq!(transcript.get_f64("gnomad_mt_heteroplasmic"), Ok(0.0123));
}

#[test]
fn test_parse_vep_freq_gnomad_vep107() {
    // GIVEN VEP 107 style gnomADg fields (trimmed from scout's vep_107_csq)
    let transcripts = transcripts(
        "Allele|Consequence|Feature|gnomADg|gnomADg_AF_POPMAX|gnomADg_AF",
        "TA|intron_variant|ENST00000429747|3:131697231-131697243|7.69231e-02|5.47945e-02",
    );

    // THEN gnomADg_AF is the gnomAD AF, and the POPMAX the gnomAD max
    assert_eq!(transcripts[0].get_f64("gnomad_maf"), Ok(0.0547945));
    assert_eq!(transcripts[0].get_f64("gnomad_max"), Ok(0.0769231));
}

#[test]
fn test_parse_vep_freq_exac() {
    // GIVEN a transcript with EXAC_MAX_AF
    let transcripts = transcripts(
        "Allele|Consequence|Feature|EXAC_MAX_AF",
        "C|missense_variant|ENST1|0.01",
    );

    // THEN both exac_max and exac_maf are set
    assert_eq!(transcripts[0].get_f64("exac_max"), Ok(0.01));
    assert_eq!(transcripts[0].get_f64("exac_maf"), Ok(0.01));
}

#[test]
fn test_parse_vep_freq_thousand_g_max() {
    // GIVEN a transcript with 1000G AFR_AF and AMR_AF frequencies
    let transcripts = transcripts(
        "Allele|Consequence|Feature|AFR_AF|AMR_AF",
        "C|missense_variant|ENST1|0.01|0.001",
    );

    // THEN thousandg_max is the highest
    assert_eq!(transcripts[0].get_f64("thousandg_max"), Ok(0.01));
}

#[test]
fn test_parse_vep_freq_skips_invalid_values() {
    // GIVEN frequency fields that are empty, ".", alphabetic, or not a number
    let transcripts = transcripts(
        "Allele|Consequence|Feature|AF|gnomAD_AF|AFR_AF|gnomAD_NFE_AF|MAX_AF_POPS",
        "C|missense_variant|ENST1|.|gnomAD_EAS||1e-3x|",
    );

    // THEN no frequencies are set
    for key in [
        "thousand_g_maf",
        "gnomad_maf",
        "thousandg_max",
        "gnomad_max",
    ] {
        assert!(!transcripts[0].contains_key(key), "{key} should not be set");
    }
}

#[test]
fn test_parse_superdups_fractmatch() {
    // GIVEN a transcript with UCSC superdups fracMatch values
    let transcripts = transcripts(
        "Allele|Consequence|Feature|genomic_superdups_frac_match",
        "C|missense_variant|ENST1|0.992904&0.98967",
    );

    // THEN all values are parsed
    assert_eq!(
        transcripts[0].get("superdups_fracmatch"),
        Some(&Bson::Array(vec![
            Bson::Double(0.992904),
            Bson::Double(0.98967)
        ]))
    );
}

#[test]
fn test_parse_cadd() {
    // GIVEN the scout vep_csq fixture with CADD_PHRED
    let transcripts = transcripts(SCN5A_CSQ_FORMAT, SCN5A_CSQ);

    // THEN the CADD score is parsed as a float
    assert_eq!(transcripts[0].get_f64("cadd"), Ok(1.179));
}

#[test]
fn test_parse_hg38_mane_transcripts() {
    // GIVEN the scout vep_csq fixture with MANE_SELECT and MANE_PLUS_CLINICAL
    let transcripts = transcripts(SCN5A_CSQ_FORMAT, SCN5A_CSQ);
    assert_eq!(transcripts.len(), 1);

    // THEN both MANE annotations are parsed
    for transcript in &transcripts {
        assert_eq!(
            transcript.get_str("mane_select_transcript"),
            Ok("NM_000335.5")
        );
        assert_eq!(
            transcript.get_str("mane_plus_clinical_transcript"),
            Ok("NM_001099404.2")
        );
    }
}

#[test]
fn test_parse_legacy_mane_field() {
    // GIVEN an older VEP annotation with a single MANE field
    let transcripts = transcripts("Allele|Feature|MANE", "C|ENST1|NM_000335.5");

    // THEN it is used as the MANE Select transcript
    assert_eq!(
        transcripts[0].get_str("mane_select_transcript"),
        Ok("NM_000335.5")
    );
    assert!(!transcripts[0].contains_key("mane_plus_clinical_transcript"));
}

#[test]
fn test_parse_dbsnp_cosmic_and_clnsig() {
    // GIVEN the scout vep_csq fixture with rs, COSV and other ids in Existing_variation
    let (transcripts, variant) = transcripts_with_variant(SCN5A_CSQ_FORMAT, SCN5A_CSQ, &[], "");

    // THEN rs ids are dbSNP and COSM/COSV ids are COSMIC, at transcript level
    assert_eq!(transcripts[0].get("dbsnp"), Some(&strings(&["rs7430407"])));
    assert_eq!(
        transcripts[0].get("cosmic"),
        Some(&strings(&["COSV61125864"]))
    );
    assert_eq!(
        transcripts[0].get("clnsig"),
        Some(&strings(&["likely_benign", "benign"]))
    );

    // AND they are collected on the variant
    assert_eq!(variant.get_str("dbsnp_id"), Ok("rs7430407"));
    assert_eq!(variant.get("cosmic_ids"), Some(&strings(&["COSV61125864"])));
}

#[test]
fn test_parse_cosmic_from_info() {
    // GIVEN COSMIC ids both in CSQ and in a separate INFO field
    let (_, variant) = transcripts_with_variant(
        "Allele|Feature|Existing_variation",
        "C|ENST1|COSV1",
        &["##INFO=<ID=COSMIC,Number=1,Type=String,Description=\"COSMIC ids\">"],
        "COSMIC=COSV1&COSV2&",
    );

    // THEN they are merged and deduplicated on the variant
    let mut cosmic_ids: Vec<&str> = variant
        .get_array("cosmic_ids")
        .unwrap()
        .iter()
        .map(|id| id.as_str().unwrap())
        .collect();
    cosmic_ids.sort();
    assert_eq!(cosmic_ids, ["COSV1", "COSV2"]);
}

#[test]
fn test_parse_existing_dbsnp_id_is_kept() {
    // GIVEN a variant that already has a dbsnp_id (e.g. from the ID column)
    let (record, field_indices) = vep_record(
        "Allele|Feature|Existing_variation",
        "C|ENST1|rs123",
        &[],
        "",
    );
    let mut variant = doc! { "dbsnp_id": "rs999" };

    // WHEN parsing the transcripts
    parse_vep_transcripts(&record, &field_indices, &mut variant);

    // THEN the existing id is not overwritten
    assert_eq!(variant.get_str("dbsnp_id"), Ok("rs999"));
}

#[test]
fn test_parse_transcripts_without_feature_are_skipped() {
    // GIVEN a regulatory annotation without a Feature id next to a transcript
    let transcripts = transcripts(
        "Allele|Consequence|Feature",
        "C|regulatory_region_variant|,C|missense_variant|ENST1",
    );

    // THEN only the transcript with an id is kept
    assert_eq!(transcripts.len(), 1);
    assert_eq!(transcripts[0].get_str("transcript_id"), Ok("ENST1"));
}

#[test]
fn test_parse_transcripts_without_csq_header() {
    // GIVEN a record without a CSQ header
    let record = crate::vcf_record::record(
        &["##contig=<ID=1,length=249250621>"],
        "1\t100\t.\tA\tG\t.\tPASS\t.",
    );
    let mut variant = Document::new();

    // THEN no transcripts are parsed
    let transcripts =
        parse_vep_transcripts(&record, &std::collections::HashMap::new(), &mut variant);
    assert!(transcripts.is_empty());
    assert!(variant.is_empty());
}
