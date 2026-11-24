# AGENTS.md

Guidance for AI coding agents (and humans) working in `scout_loader`.

## What this project is

`scout_loader` is a Rust CLI that parses VCF files for a case and loads the
variants into a [Scout](https://github.com/Clinical-Genomics/scout) MongoDB
database. It is a fast replacement for the variant *parse/build/load* path of
Scout (`scout/parse/variant/`, `scout/build/variant/`). It has no web UI;
Scout itself handles display. Scout downloads and runs the released binaries
(see `RELEASE.md`).

**The documents it writes must match the ones Scout's Python loader would
write**: same keys, same types, same "missing value" behaviour. Scout reads
them back, so a renamed key, or an `Int64` where Scout writes an `Int32`, is a
bug even if every test here passes.

## Porting features from Scout

Most work is porting a Scout PR. For each one:

1. Read the Scout diff (`gh pr diff <N> -R Clinical-Genomics/scout`). A local
   checkout may be at `../scout`.
2. Only port the parse/build parts. Ignore templates, controllers, and the
   front end.
3. Scout PRs often bundle refactors. Check whether `scout_loader` already
   matches before porting them; it often does.
4. Copy Scout's output contract exactly: the document key names, which values
   count as missing (Scout often leaves the key out rather than writing
   `null`), and the BSON type.
5. Port or adapt Scout's tests. Mention the Scout PR in the loader PR.

Rough mapping from Scout modules to loader modules:

| Scout (Python)                         | scout_loader (Rust)                                     |
| -------------------------------------- | ------------------------------------------------------- |
| `scout/parse/variant/variant.py`       | `src/parse/vcf.rs` (`process_vcf`) + `src/parse/*.rs`   |
| `scout/parse/variant/genotype.py`      | `src/parse/genotypes.rs`                                |
| `scout/build/variant/genotype.py`      | the `gt_obj` / `format_tag` list in `parse_genotype`    |
| `scout/parse/variant/transcript.py`, VEP CSQ | `src/parse/vep/`                                  |
| `scout/build/variant/gene.py`, `transcript.py` | `src/build/gene.rs`, `src/build/transcript.rs`  |
| `scout/parse/variant/frequency.py`     | `src/parse/frequencies.rs`, `src/parse/vep/frequencies.rs` |
| `scout/parse/variant/compound.py`      | `src/parse/compounds.rs` (incl. `update_compounds`)     |
| `scout/parse/variant/models.py`        | `src/parse/genetic_models.rs`                           |
| `scout/parse/variant/clnsig.py`        | `src/parse/vep/clnsig.rs`, `src/parse/onco_clnsig.rs`   |

## Layout

- `src/main.rs`: CLI (clap). Reads the case YAML and calls `parser::parse`.
- `src/parser.rs`: picks VCFs by category, clinical or research.
- `src/parse/vcf.rs`: `process_vcf`, the per-record loop that builds the
  variant `Document`. Category-specific annotations go in the `match category`
  block.
- `src/parse/`: one module per annotation area. Most are small functions of
  the form `fn(record: &Record, …) -> Option<T>`.
- `src/loader.rs`: MongoDB access (genes, panels, causatives, bulk inserts).
- `src/updater.rs`: post-load variant rank updates (`update_variant_rank`).
  Compound updates are in `src/parse/compounds.rs`.
- `src/models/`: case config, sample, and variant category types.
- `tests/`: integration-style tests. `tests/parse.rs` pulls in
  `tests/parse/*.rs` through `#[path]` modules. **Add a new test file there**
  or it will not be compiled.

## Conventions

- Read FORMAT and INFO values with `rust_htslib`: `record.format(b"TAG").integer()`
  or `.float()` / `.string()`, and `record.info(b"TAG")`. Missing integers come
  back as negative sentinels (`i32::MIN`, vector end), so filter with `>= 0`
  where a value cannot be negative.
- If an optional field is absent, don't insert its key. Use
  `if let Some(v) = … { doc.insert(…) }`, which matches Scout's behaviour.
- Genotype-level FORMAT fields: parse in `parse_genotype`, insert into
  `gt_call`, then add the key to the `format_tag` allow-list so it is copied
  into `gt_obj`. Keep that list in sync with Scout's `build_genotype`.
- Reuse the helpers in `genotypes.rs` (`parse_format_entry_single_integer`,
  `_float`, `_string`, `parse_format_entry`) and in `info.rs`
  (`insert_info_*`, `parse_info_*`) before writing new ones.
- Add doc comments to public parse functions. Include the VCF header line
  (`##FORMAT=<ID=…>`) for any new tag you handle.
- Put tests under `tests/`, not in `#[cfg(test)]` modules in `src/`, so they
  are easy to find. Make a function `pub` (with a doc comment) when a test
  needs it. Keep test-only constants in the test file, at the top after the
  imports.
- Tests use given/when/then comments. Build records from inline VCF text with
  the shared fixtures in `tests/fixtures/vcf_record.rs` (`record`,
  `records_from_vcf`, `vep_record`) instead of writing temp files in each
  test; see `trgt_record` in `tests/parse/genotypes.rs`. htslib prints
  "Contig not defined" warnings for header-less test VCFs. These are harmless.
- Rust edition 2024. Let-chains (`if let … && let …`) are used throughout.

## Commands

```bash
cargo build
cargo test --locked --all-targets           # some tests need MongoDB on 127.0.0.1:27017
cargo test --test parse genotypes           # one test module
cargo fmt --all --check                     # CI lint
cargo clippy --locked --all-targets -- -D warnings   # CI lint, warnings are errors
```

Some tests in `tests/loader.rs` and `tests/mongodb_smoke.rs` need a local
MongoDB, for example `docker run -d -p 27017:27017 mongo:7`. CLI tests set
`TEST_ENV=1` to skip the institute lookup.

## Before opening a PR

- `cargo fmt`, clippy with `-D warnings`, and the full test suite all pass.
- Add an entry under `## [unreleased]` in `CHANGELOG.md` (Added / Changed /
  Fixed) with the PR number. CI checks for this.
- Don't bump the version in `Cargo.toml`. That happens at release time
  (`RELEASE.md`).
