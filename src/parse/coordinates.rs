use rust_htslib::bcf::header::HeaderView;
use rust_htslib::bcf::Record;

pub fn parse_coordinates(
    record: &Record,
    header: &HeaderView,
) -> (String, u64, u64) {
    let rid = record.rid().expect("missing chromosome");

    let chromosome = String::from_utf8_lossy(
        header.rid2name(rid).expect("unknown chromosome")
    )
    .to_string();

    let position = (record.pos() + 1) as u64;

    let end = position;

    (chromosome, position, end)
}