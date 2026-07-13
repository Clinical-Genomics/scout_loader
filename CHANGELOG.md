# Change Log
All notable changes to this project will be documented in this file.
This project adheres to [Semantic Versioning](http://semver.org/).

About changelog [here](https://keepachangelog.com/en/1.0.0/)

## [unreleased]
### Fixed
- Normalize chromosome names by stripping the chr prefix (#4)
### Added
- Use the `rust-htslib` for reading VCF files (#1)
- First Variant model and parsing modules (#2)
- Parse alleles and missing docstrings (#3)
- Variant type enum (clinical or research) key/values passed to main function (#5)
- Custom issues and pull requests templates (#9)
- Pass case _id from CLI (#10)
- Parse variant.FILTER field (#11)
- Parse variant.QUAL field (#12)