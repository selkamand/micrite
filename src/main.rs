use clap::{Parser, Subcommand};
use env_logger::Env;
use std::fs::read_to_string;
use std::path::PathBuf;

#[derive(Subcommand)]
enum Commands {
    /// Screen bam for microbial presense
    Screen {
        /// File with paths to bam files (newline separated)
        #[arg(short, long, value_name = "BAM File")]
        bam: PathBuf,

        /// Path to Kraken Database
        #[arg(long, value_name = "Kraken Database")]
        db_kraken: PathBuf,

        /// Threads
        #[arg(short, long, default_value_t = 1)]
        threads: u8,

        /// Confidence Threshold
        #[arg(short, long, default_value_t = 0.5)]
        confidence: f32,

        /// Output std file
        #[arg(long, default_value_t = false)]
        cleanup_std_file: bool,

        /// Delete unmapped reads extracted from bam file after use
        #[arg(long, default_value_t = false)]
        cleanup_unmapped: bool,

        /// Include zero counts in kraken report
        #[arg(long, default_value_t = false)]
        report_zero_counts: bool,

        /// Minimum proportion of unmapped reads that must be classified as a microbe to flag as a hit
        #[arg(short = 'p', long, default_value_t = 0.01)]
        min_prop_unmapped_reads: f32,

        /// Minimum number of unmapped reads that must be classified as a microbe to flag as a hit
        #[arg(short = 'n', long, default_value_t = 50)]
        min_number_unmapped_reads: u64,

        /// Only flag presense of known oncogenic viruses
        #[arg(short = 'O', long, default_value_t = false)]
        oncogenic_only: bool,

        // DEACON settings
        /// Path to deacon host database used for host-depletion.
        /// We advise use of the 3.4gb precompiled panhuman-1 index available from https://github.com/bede/deacon.
        #[arg(long, value_name = "Deacon Database")]
        db_host: PathBuf,

        /// Minimum absolute number of minimizer hits for a match
        #[arg(short, long, default_value_t = 2)]
        absolute_threshold: u8,

        /// Minimum relative proportion (0.0-1.0) of minimizer hits for a match
        #[arg(short, long, default_value_t = 0.01)]
        relative_threshold: f32,
    },

    /// Validate reads truly come from a microbe of interest
    Sleuth,

    /// Subtype a genome
    Subtype,
}

#[derive(Parser)]
#[command(version, about, long_about = None)]
struct Cli {
    #[command(subcommand)]
    command: Commands,

    /// Output directory
    #[arg(short, long, value_name = "OUTDIR")]
    outdir: PathBuf,
}

fn setup_logger() {
    env_logger::Builder::from_env(Env::default().default_filter_or("info"))
        .format_timestamp_secs()
        .init();
}
fn main() {
    // Setup logger
    setup_logger();

    // Run main function
    if let Err(err) = run() {
        // One place to log the top-level error (and its causes)
        log::error!("{err:?}"); // or `{:#}` for pretty, or just `{}` for brief
        std::process::exit(1);
    }
}

fn run() -> Result<(), anyhow::Error> {
    let cli = Cli::parse();

    match &cli.command {
        Commands::Screen {
            bam,
            db_host,
            absolute_threshold,
            relative_threshold,
            db_kraken,
            threads,
            confidence,
            cleanup_std_file,
            cleanup_unmapped,
            report_zero_counts,
            min_prop_unmapped_reads,
            min_number_unmapped_reads,
            oncogenic_only,
        } => {
            let kraken_config = micrite::kraken::KrakenConfig {
                krakendb: db_kraken.clone(),
                threads: *threads,
                confidence: confidence.to_string(),
                cleanup_std_file: *cleanup_std_file,
                cleanup_unmapped: *cleanup_unmapped,
                report_zero_counts: *report_zero_counts,
                kraken_hit_thresholds: micrite::kraken::KrakenHitThresholds {
                    min_prop_unmapped_reads: *min_prop_unmapped_reads,
                    min_number_reads: *min_number_unmapped_reads,
                    oncogenic_only: *oncogenic_only,
                },
                outdir: cli.outdir.display().to_string(),
            };
            let deacon_config = micrite::hostdepletion::DeaconConfig {
                db: db_host.clone(),
                relative_threshold: *relative_threshold,
                absolute_threshold: *absolute_threshold,
            };

            // Identify Microbes from BAM
            micrite::bam::bam2microbes(bam, &kraken_config, &deacon_config)?;
        }

        Commands::Sleuth => panic!("Validation is not yet implemented"),
        Commands::Subtype => todo!("Subtyping is not yet implemented"),
    }

    log::info!("Micrite run completed");
    Ok(())
}

fn _read_lines(filename: &str) -> Vec<String> {
    read_to_string(filename)
        .unwrap() // panic on possible file-reading errors
        .lines() // split the string into an iterator of string slices
        .map(String::from) // make each slice into a string
        .collect() // gather them together into a vector
}
