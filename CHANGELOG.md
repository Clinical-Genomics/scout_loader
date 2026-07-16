# Change Log
All notable changes to this project will be documented in this file.
This project adheres to [Semantic Versioning](http://semver.org/).

About changelog [here](https://keepachangelog.com/en/1.0.0/)

## [unreleased]
### Added
- Use the `rust-htslib` for reading VCF files (#1)
- First Variant model and parsing modules (#2)
- Parse alleles and missing docstrings (#3)
- Variant type enum (clinical or research) key/values passed to main function (#5)
- Custom issues and pull requests templates (#9)
- Pass case _id from CLI (#10)
- Parse variant.FILTER field (#11)
- Parse variant.QUAL field (#12)
- Parse compounds field (#14)
- Parse the 4 variant IDs: `simple_id`, `variant_id`, `display_name`, `document_id` (#15)
- Parse `rank_score` and `norm_rank_score`from VCF (#16)
- Infer complete coordinate data (#17)
- Parse cytobands and collect cytoband start and end for each variant (#18)
- Parse genmod genetic models (#19)
### Fixed
- Normalize chromosome names by stripping the chr prefix (#4)