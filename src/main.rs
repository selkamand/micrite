use clap::{Parser, Subcommand};
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
        #[arg(short, long, value_name = "Kraken Database")]
        db: PathBuf,

        /// Threads
        #[arg(short, long, default_value_t = 1)]
        threads: u8,

        /// Confidence Threshold
        #[arg(short, long, default_value_t = 0.5)]
        confidence: f32,

        /// Output std file
        #[arg(long, default_value_t = false)]
        output_std_file: bool,

        // Minimum proportion of unmapped reads that must be classified as a microbe to flag as a hit
        #[arg(short = 'p', long, default_value_t = 0.01)]
        min_prop_unmapped_reads: f32,

        // Minimum number of unmapped reads that must be classified as a microbe to flag as a hit
        #[arg(short = 'n', long, default_value_t = 50)]
        min_number_unmapped_reads: u64,

        /// Only flag presense of known oncogenic viruses
        #[arg(short = 'O', long, default_value_t = false)]
        oncogenic_only: bool,
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

fn main() {
    let cli = Cli::parse();

    match &cli.command {
        Commands::Screen {
            bam,
            db,
            threads,
            confidence,
            output_std_file,
            min_prop_unmapped_reads,
            min_number_unmapped_reads,
            oncogenic_only,
        } => {
            let config = micrite::kraken::KrakenConfig {
                krakendb: db.clone(),
                threads: *threads,
                confidence: confidence.to_string(),
                output_std_file: *output_std_file,
                kraken_hit_thresholds: micrite::kraken::KrakenHitThresholds {
                    min_prop_unmapped_reads: *min_prop_unmapped_reads,
                    min_number_reads: *min_number_unmapped_reads,
                    oncogenic_only: *oncogenic_only,
                },
                outdir: cli.outdir.display().to_string(),
            };

            // Identify Microbes from BAM
            micrite::bam::bam2microbes(bam, &config);
        }
        Commands::Sleuth => panic!("Validation is not yet implemented"),
        Commands::Subtype => todo!("Subtyping is not yet implemented"),
    }

    eprintln!("finished")
}

fn _read_lines(filename: &str) -> Vec<String> {
    read_to_string(filename)
        .unwrap() // panic on possible file-reading errors
        .lines() // split the string into an iterator of string slices
        .map(String::from) // make each slice into a string
        .collect() // gather them together into a vector
}
