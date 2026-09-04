# Scout Loader

[![Tests](https://github.com/Clinical-Genomics/scout_loader/actions/workflows/tests_n_coverage.yml/badge.svg)](https://github.com/Clinical-Genomics/scout_loader/actions/workflows/tests_n_coverage.yml)
[![codecov](https://codecov.io/gh/Clinical-Genomics/scout_loader/graph/badge.svg)](https://codecov.io/gh/Clinical-Genomics/scout_loader)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

`scout_loader` is a command-line tool intended to load **clinical and research variants from VCF files into a Scout database**.

The repository is designed to work alongside the [Scout](https://github.com/Clinical-Genomics/scout) repository and uses the case configuration to determine which VCF files are available for a case.

## Requirements

* A running and accessible **Scout MongoDB database**
* A Scout case configuration YAML file containing the VCF files and case information
* A MongoDB configuration file

The MongoDB configuration can be provided explicitly using `--config`. If no configuration file is specified, `scout_loader` looks for the default `config.toml` in the current directory.

## MongoDB configuration

`scout_loader` requires access to a running Scout MongoDB database. The database connection is configured using a TOML file containing the MongoDB URI and database name.

A minimal configuration looks like:

```toml
mongo_uri = "mongodb://127.0.0.1:27017"
mongo_dbname = "scout-demo"
```

* `mongo_uri` — MongoDB connection URI.
* `mongo_dbname` — Name of the Scout database to load variants into.

By default, `scout_loader` looks for a file named `config.toml` in the current working directory.

Alternatively, a different configuration file can be provided using the `--config` option:

```bash
scout_loader \
    --config /path/to/config.toml \
    --case-config case.yaml
```

For example, the repository's demo configuration can be used to connect to a local `scout-demo` database:

```toml
mongo_uri = "mongodb://127.0.0.1:27017"
mongo_dbname = "scout-demo"
```

The MongoDB database must be available before running `scout_loader`.

## Command-line usage

The basic command is:

```bash
scout_loader --case-config <CASE_CONFIG>
```

For example:

```bash
scout_loader \
    --case-config case.yaml
```

If the MongoDB configuration is not in the default `config.toml` location, provide it explicitly:

```bash
scout_loader \
    --config /path/to/config.toml \
    --case-config case.yaml
```

### Command-line options

| Option                      | Required | Description                                                             |
| --------------------------- | -------- | ----------------------------------------------------------------------- |
| `--case-config <PATH>`      | Yes      | Path to the case configuration YAML file.                               |
| `--config <PATH>`           | No       | Path to the MongoDB configuration TOML file. Defaults to `config.toml`. |
| `--categories <CATEGORIES>` | No       | Load only the specified variant categories.                             |
| `--research`                | No       | Load research VCFs instead of clinical VCFs.                            |

## VCF selection

The case configuration specifies the VCF files available for the case. VCFs are organized by **variant category**, for example:

* `snv`
* `sv`
* `str`
* `mei`
* `cancer`
* `cancer_sv`
* `fusion`

Both clinical and research VCFs can be specified independently. For example:

```yaml
vcf_snv: /path/to/case.clinical.vcf.gz
vcf_sv: /path/to/case.clinical.SV.vcf.gz
vcf_snv_research: /path/to/case.research.vcf.gz
vcf_sv_research: /path/to/case.research.SV.vcf.gz
```

### Default behavior

When only `--case-config` is provided, `scout_loader` loads **all available clinical VCFs** from the configuration.

The normal clinical filtering and variant-loading rules are applied during parsing.

For example:

```bash
scout_loader --case-config case.yaml
```

loads the available clinical VCFs.

### Loading research VCFs

Use `--research` to load **research VCFs instead of clinical VCFs**:

```bash
scout_loader \
    --case-config case.yaml \
    --research
```

When `--research` is specified, only research VCFs are selected.

### Selecting specific categories

The `--categories` option can be used to restrict loading to specific **variant categories**.

For example, to load only SNVs and SVs:

```bash
scout_loader \
    --case-config case.yaml \
    --categories snv,sv
```

The same category selection can be combined with `--research`:

```bash
scout_loader \
    --case-config case.yaml \
    --categories snv,sv \
    --research
```

This loads only the **research SNV and research SV VCFs**.

In other words:

| Command                                                               | VCFs loaded                 |
| --------------------------------------------------------------------- | --------------------------- |
| `scout_loader --case-config case.yaml`                                | All available clinical VCFs |
| `scout_loader --case-config case.yaml --research`                     | All available research VCFs |
| `scout_loader --case-config case.yaml --categories snv,sv`            | Clinical SNV and SV VCFs    |
| `scout_loader --case-config case.yaml --categories snv,sv --research` | Research SNV and SV VCFs    |

Only VCFs that are actually present in the case configuration are loaded.

## Case configuration

The case configuration YAML provides the information required to load the case, including the case identifier, samples, genome build, and available VCFs.

For example:

```yaml
owner: test_institute
family: case_123

human_genome_build: "37"

samples:
  - sample_id: SAMPLE1
    sample_name: sample-one

vcf_snv: /path/to/case.clinical.vcf.gz
vcf_sv: /path/to/case.clinical.SV.vcf.gz
vcf_snv_research: /path/to/case.research.vcf.gz
vcf_sv_research: /path/to/case.research.SV.vcf.gz
```

Additional Scout configuration fields may be present in the YAML file; `scout_loader` only reads the fields required for loading the case and variants.

