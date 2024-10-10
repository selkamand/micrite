use std::{path::PathBuf, string};

pub struct KrakenConfig {
    pub krakendb: PathBuf,
    pub threads: u8,
    pub confidence: String,
    pub outdir: String,
}
pub fn run_kraken(fasta: std::path::PathBuf, config: KrakenConfig) {
    std::fs::create_dir_all(&config.outdir).expect("Failed to create output directory");
    let filename = fasta.file_stem().expect("Failed to extract fasta file stem (are you sure you supplied a filepath and not a directory?)").to_str().expect("failed filepath to str conversion");
    let outfile_prefix = format!("{}/{}", config.outdir, filename);
    let outfile_report = format!("{}.kreport", outfile_prefix);
    let outfile_unclassified = format!("{}.unclassified", outfile_prefix);
    let outfile_classified = format!("{}.classified", outfile_prefix);
    let outfile_output = format!("{}.output.tsv", outfile_prefix);
    let kraken_command = which::which("kraken2")
        .expect("Kraken2 not found. Please ensure it is installed and added to your PATH.");

    let db: std::borrow::Cow<'_, str> =
        shellexpand::full(config.krakendb.to_str().expect("failed to_str()"))
            .expect("Failed expansion of DB filepath");

    eprintln!("Running Kraken");
    let output = std::process::Command::new(kraken_command)
        .args(["--db", db.as_ref()])
        .args(["--threads", &config.threads.to_string()])
        .args(["--confidence", &config.confidence])
        .args(["--unclassified-out", &outfile_unclassified])
        .args(["--classified-out", &outfile_classified])
        .args(["--output", &outfile_output])
        .args(["--report", &outfile_report])
        .arg(fasta)
        .output()
        .expect("Failed to run Kraken2 classification");

    if !output.status.success() {
        let stderr_str = String::from_utf8_lossy(&output.stderr);
        panic!(
            "\tKraken Run Failed. Stderr\n========\n{}\n========",
            stderr_str
        )
    }
    // .arg(arg);
}
