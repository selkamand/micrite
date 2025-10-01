use anyhow::Context;
use std::path::{Path, PathBuf};

/// Configuration for running Deacon host read depletion.
pub struct DeaconConfig {
    /// Path to the Deacon minimiser index (e.g., panhuman-1).
    pub db: PathBuf,
    /// `-r/--rel-threshold`: minimum relative proportion (0.0–1.0) of minimizer hits for a match.
    /// Typical default is `0.01` (1%).
    pub relative_threshold: f32,
    /// `-a/--abs-threshold`: minimum absolute number of minimizer hits for a match.
    /// Typical default is `2`.
    pub absolute_threshold: u8,
}

/// Run host read depletion via [`deacon`](https://github.com/bede/deacon) using its
/// `filter` subcommand with **invert filtering** (`-d/--deplete`), i.e. **discard host-matching
/// sequences** and keep the non-host reads.
///
/// Internally, this executes an equivalent of:
/// ```text
/// deacon filter -d -a <ABS_THRESHOLD> -r <REL_THRESHOLD> -o <OUTPUT> <DB> <FASTA>
/// ```
///
/// # Parameters
///
/// - `fasta`: Path to the input FASTA/FASTQ file to deplete (single-end).
/// - `fasta_output`: Destination path for the **non-host** output FASTA/FASTQ. Compression is
///   auto-detected by extension (e.g., `.gz`, `.zst`) by Deacon.
/// - `config`: Thresholds and database path (see [`DeaconConfig`]).
///
/// # Behavior
///
/// - Uses `-d/--deplete` so **matching (host) reads are removed**.
/// - Forwards:
///   - `-a/--abs-threshold` from `config.absolute_threshold`
///   - `-r/--rel-threshold` from `config.relative_threshold`
///   - `-o/--output` to `fasta_output`
///   - `<DB>` from `config.db` and `<FASTA>` from `fasta`
/// - Captures child stdout/stderr; stdout (if any) is logged at `debug`, and non-zero exit codes
///   include stderr in the error message.
/// - This wrapper targets single-end data. Paired-end output (`-O/--output2`) isn’t wired here.
///
/// # Returns
///
/// On success, returns `fasta_output` as a `PathBuf`.
///
/// # Errors
///
/// Returns an error if:
/// - `deacon` isn’t installed or not on `PATH`.
/// - `config.db` does not exist.
/// - The process fails to spawn or exits with a non-zero status (stderr included).
///
/// # Notes
///
/// - Common Deacon defaults (if you want to mirror them in your config):
///   - `relative_threshold = 0.01` (1%)  
///   - `absolute_threshold = 2`
/// - Additional useful Deacon flags not exposed here:  
///   `-t/--threads`, `-s/--summary`, `-p/--prefix-length`, `--compression-level`,
///   `-q/--quiet`, `--debug`, `-R/--rename`.
/// - Ensure `relative_threshold` is within `0.0..=1.0` for valid runs.
///
/// ```
pub fn host_depletion(
    fasta: &Path,
    fasta_output: &Path,
    config: &DeaconConfig,
) -> Result<PathBuf, anyhow::Error> {
    // Locate `deacon` in PATH
    let deacon_command = which::which("deacon")
        .context("`deacon` not found. Ensure it is installed and in your PATH. See https://github.com/bede/deacon")?;

    // Verify index exists
    if !&config.db.exists() {
        anyhow::bail!(
            "Failed to find deacon minimiser index: {}",
            &config.db.display()
        );
    }

    // Get Threshold info
    let a = config.absolute_threshold.to_string();
    let r = config.relative_threshold.to_string();

    // Build command
    let mut cmd = std::process::Command::new(deacon_command);
    cmd.arg("filter")
        .arg("-d")
        .args(["-a", &a])
        .args(["-r", &r])
        .args([
            "-o",
            fasta_output
                .to_str()
                .context("Failed to convert fasta_output to str")?,
        ])
        .arg(config.db.clone())
        .arg(fasta);

    log::info!("Running Deacon: {cmd:?}");

    // Run and capture output
    let output = cmd
        .output()
        .context("Failed to run deacon host depletion")?;

    // Log stdout (useful for diagnostics at higher verbosity)
    if !output.stdout.is_empty() {
        let stdout_str = String::from_utf8_lossy(&output.stdout);
        log::debug!("Deacon stdout:\n{stdout_str}");
    }

    // Handle failures with helpful stderr
    if !output.status.success() {
        let stderr_str = String::from_utf8_lossy(&output.stderr);
        anyhow::bail!(
            "Deacon run failed (exit code: {}).\n--- STDERR ---\n{stderr_str}\n---------------",
            output.status.code().unwrap_or(-1)
        );
    }

    log::info!(
        "Deacon non-host reads written to {}",
        fasta_output.display()
    );
    Ok(fasta_output.to_path_buf())
}
