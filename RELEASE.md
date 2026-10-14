# Release procedure

This document describes how to prepare and publish a new release of `scout_loader`.

## 1. Create a release branch

Create a new branch from the appropriate base branch, using the release version in the branch name.

For example, for version `0.1.2`:

```bash
git checkout <base-branch>
git pull
git checkout -b version_0.1.2
```

## 2. Update the changelog

Open `CHANGELOG.md` and update the `[unreleased]` section to reflect the new release version.

Replace the `[unreleased]` heading with the new version number and release date, following the existing changelog format.

For example:

```markdown
## [0.1.2]
```

Ensure that all changes intended for this release are documented in the corresponding section.

## 3. Update dependencies

Run the following command to update the Rust dependencies recorded in `Cargo.lock`:

```bash
cargo update
```

Run the relevant tests to verify that the project still builds and works correctly.

## 4. Update the version in the TOML file

Update the package version in `Cargo.toml` to match the intended release version.

For example:

```toml
[package]
version = "0.1.2"
```

Make sure the version is consistent across `Cargo.toml`, `CHANGELOG.md`, and the planned GitHub release tag.

## 5. Commit and push the release branch

Review all modified files, including the changelog, dependency lockfile, and TOML configuration.

Stage and commit the changes, then push the release branch to GitHub:

```bash
git status
git add CHANGELOG.md Cargo.toml Cargo.lock
git commit -m "Prepare release v0.1.2"
git push
```

If additional files were modified as part of the release preparation, include them in the commit as appropriate.

Open a pull request for the release branch and merge it following the project's normal review process.

## 6. Create the GitHub release

After the release branch has been merged, create the release manually on GitHub.

1. Open the [scout_loader GitHub Releases page](https://github.com/Clinical-Genomics/scout_loader/releases).
2. Create a new release using the appropriate target branch or commit.
3. Create a tag matching the release version, prefixed with `v`.

   For example: `v0.1.2`

4. Add the release title and notes, following the project's existing release conventions.
5. Publish the release to trigger the associated GitHub Actions release workflow.

## 7. Verify the release workflow and assets

**Do not consider the release complete until the GitHub Actions workflow has finished successfully and the expected assets are available.**

After publishing the release:

1. Open the repository's **Actions** tab and locate the workflow associated with the release process.
2. Verify that the workflow was triggered for the correct tag and completed successfully.
3. Open the published release and inspect its assets.
4. Confirm that all expected release files are attached, especially the compiled binary files required to download and execute `scout_loader`.
5. Check that the expected binaries for the supported platforms and architectures are present and named correctly.
6. If any expected files are missing, investigate the workflow and its build or upload steps before considering the release ready for use.

The release assets are particularly important because Scout downloads and executes the appropriate `scout_loader` binary. A published GitHub release without the expected executable files can prevent the software from running correctly.

## Release checklist

- [ ] Release branch created with the correct version.
- [ ] `CHANGELOG.md` updated.
- [ ] Dependencies updated using `cargo update`.
- [ ] Version updated in `Cargo.toml`.
- [ ] Modified files reviewed, committed, and pushed.
- [ ] Pull request reviewed and merged.
- [ ] GitHub release published with the correct version tag.
- [ ] Release GitHub Actions workflow completed successfully.
- [ ] All expected release assets, especially executable binaries, are available and ready to download.