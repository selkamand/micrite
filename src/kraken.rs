use anyhow::Context;
use std::path::PathBuf;

pub struct KrakenConfig {
    pub krakendb: PathBuf,
    pub threads: u8,
    pub confidence: String,
    pub cleanup_std_file: bool, // Should std kraken tsv mapping readnames to taxids be output (large files, but required for pulling out taxid-specific reads)
    pub cleanup_unmapped: bool, // Should unmapped reads extracted from bams be kept after use?
    pub report_zero_counts: bool, // Should kraken report include species with no read support?
    pub kraken_hit_thresholds: KrakenHitThresholds,
    pub outdir: String,
}

pub struct Microbe {
    pub name: String,
    pub taxid: String,
}

pub struct CancerMicrobes {
    microbes: Vec<Microbe>,
}

impl CancerMicrobes {
    // Check if InterestingContigs contain a particular contig name
    fn contains(&self, taxid: &str) -> bool {
        let taxids_in_set: Vec<&str> = self.microbes.iter().map(|c| c.taxid.as_str()).collect();
        taxids_in_set.contains(&taxid)
    }

    // If Taxid
    fn _taxid_to_name(&self, taxid: &str) -> Option<&str> {
        let species = self
            .microbes
            .iter()
            .filter(|c| c.taxid.as_str() == taxid)
            .map(|c| c.name.as_str())
            .next();

        species
    }
}

pub fn cancer_microbes() -> CancerMicrobes {
    // Define the microbe data as tuples of (name, taxid)
    let microbe_data = vec![
        ("Human gammaherpesvirus 8", "37296"),
        ("Human gammaherpesvirus 4 (EBV)", "10376"),
        ("Human betaherpesvirus 6A", "32603"),
        ("Human betaherpesvirus 6B", "32604"),
        ("Human betaherpesvirus 7", "10372"),
        ("Primate T-lymphotropic virus 1", "194440"),
        ("Primate T-lymphotropic virus 2", "194441"),
        ("Human papillomavirus", "10566"),
        ("Hepatitis B virus", "10407"),
        ("Hepacivirus C", "11103"),
        ("Merkel cell polyomavirus", "493803"),
        ("Betapolyomavirus macacae", "1891767"),
        ("Betapolyomavirus secuhominis", "1891763"),
        ("Betapolyomavirus hominis", "1891762"),
        ("Cytolomegalovirus", "10358"),
        ("Alphatorquevirus", "687331"),
    ];

    // Convert each tuple into a Microbe struct
    let microbes = microbe_data
        .into_iter()
        .map(|(name, taxid)| Microbe {
            name: name.to_string(),
            taxid: taxid.to_string(),
        })
        .collect();

    CancerMicrobes { microbes }
}

pub struct KrakenOutputPaths {
    pub kout: Option<PathBuf>,
    pub kreport: PathBuf,
    pub input_fasta: PathBuf,
    pub prefix: String,
}

pub fn run_kraken(
    fasta: std::path::PathBuf,
    config: &KrakenConfig,
) -> Result<KrakenOutputPaths, anyhow::Error> {
    std::fs::create_dir_all(&config.outdir).context("Failed to create output directory")?;
    let filename = fasta.file_stem().context("Failed to extract fasta file stem (are you sure you supplied a filepath and not a directory?)")?.to_str().context("failed filepath to str conversion")?;
    let outfile_prefix = format!("{}/{}", config.outdir, filename);
    let outfile_report = format!("{outfile_prefix}.kreport");
    // let outfile_unclassified = format!("{}.unclassified", outfile_prefix);
    // let outfile_classified = format!("{}.classified", outfile_prefix);
    let outfile_output = match config.cleanup_std_file {
        false => format!("{outfile_prefix}.kout.tsv"),
        true => "-".to_string(),
    };

    let kraken_command = which::which("kraken2")
        .context("Kraken2 not found. Please ensure it is installed and added to your PATH.")?;

    let db: std::borrow::Cow<'_, str> =
        shellexpand::full(config.krakendb.to_str().context("failed to_str()")?)
            .context("Failed expansion of DB filepath")?;

    log::info!("\nRunning Kraken:");

    // Build KrakenCommand
    let mut binding = std::process::Command::new(kraken_command);
    let cmd_kraken = binding
        .args(["--db", db.as_ref()])
        .args(["--threads", &config.threads.to_string()])
        .args(["--confidence", &config.confidence])
        // .args(["--unclassified-out", &outfile_unclassified])
        // .args(["--classified-out", &outfile_classified])
        .args(["--output", outfile_output.as_str()])
        .args(["--report", &outfile_report])
        .arg(&fasta);

    if config.report_zero_counts {
        cmd_kraken.args(["--report-zero-counts"]);
    }
    log::info!("\nRunning Kraken: {cmd_kraken:?}");

    // Run Kraken
    let output = cmd_kraken
        .output()
        .context("Failed to run Kraken2 classification")?;

    if !output.status.success() {
        let stderr_str = String::from_utf8_lossy(&output.stderr);
        panic!("\tKraken Run Failed. Stderr\n========\n{stderr_str}\n========")
    }
    log::info!("\tKraken report saved to: {outfile_report}");

    let kout_path: Option<PathBuf> = match config.cleanup_std_file {
        true => None,
        false => Some(outfile_output.into()),
    };

    // Return the output paths
    Ok(KrakenOutputPaths {
        kout: kout_path,
        input_fasta: fasta,
        kreport: outfile_report.into(),
        prefix: outfile_prefix,
    })
}

#[derive(Debug, serde::Deserialize)]
struct KreportRecord {
    clade_percent_classified: f32, //% of reads classified as this taxid or child taxids
    clade_nreads_classified: u64,  //How many reads classified as this taxid or child taxids
    _taxon_nreads_classified: u64, //How many reads assigned directly to this taxid
    rank: String, // A rank code, indicating (U)nclassified, (R)oot, (D)omain, (K)ingdom, (P)hylum, (C)lass, (O)rder, (F)amily, (G)enus, or (S)pecies. Taxa that are not at any of these 10 ranks have a rank code that is formed by using the rank code of the closest ancestor rank with a number indicating the distance from that rank. E.g., "G2" is a rank code indicating a taxon is between genus and species and the grandparent taxon is at the genus rank.
    taxid: String, //NCBI taxonomic ID
    name: String, //Scientific Name
}

#[derive(serde::Serialize)]
struct KrakenHit<'a> {
    taxid: &'a str,
    rank: &'a str,
    name: &'a str,
    clade_percent_classified: &'a f32,
    clade_nreads_classified: &'a u64,
    oncogenic: &'a bool,
}

// struct KrakenHits {
//     hits: Vec<KrakenHits>,
//     oncogenic_only: bool, // was this oncogenics only
// }

pub struct KrakenHitThresholds {
    pub min_prop_unmapped_reads: f32,
    pub min_number_reads: u64,
    pub oncogenic_only: bool, // Only identify hits from a list of 'oncogenic' microbes. This helps reduce noise.
}

pub fn identify_kraken_hits_from_kreport(
    paths: KrakenOutputPaths,
    thresholds: &KrakenHitThresholds,
) -> Result<(), anyhow::Error> {
    // Create reader for kraken report
    let mut rdr_kreport = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .trim(csv::Trim::All)
        .from_path(paths.kreport)
        .context("failed to read kreport")?;

    // Create writer
    let oncogenic_microbe_counts: PathBuf = format!("{}.krakenhits.csv", paths.prefix).into();
    let mut wtr = csv::Writer::from_path(&oncogenic_microbe_counts)
        .context("Failed to create writer to oncogenic microbe count filepath")?;

    // List of oncogenic microbes
    let cancer_microbes = cancer_microbes();

    // Print out threshold information
    let oncogenic_only_text = match thresholds.oncogenic_only {
        true => " (oncogenic only) ",
        false => " ",
    };

    log::info!(
        "Checking kraken reports for microbes{}with >= {} supporting reads & account for >= {:4.1} % of all unmapped reads",
        oncogenic_only_text,
        thresholds.min_number_reads,
        thresholds.min_prop_unmapped_reads
    );

    let mut n_non_oncogenics_excluded: u64 = 0;
    let mut n_microbial_hits: u64 = 0;
    for records_result in rdr_kreport.deserialize() {
        let record: KreportRecord = records_result.context("Failed to read kreport record")?;

        let is_oncogenic_microbe = cancer_microbes.contains(record.taxid.as_str());

        // TODO: - add the option to normalise read counts based on a reference matrix.

        if record.clade_nreads_classified > thresholds.min_number_reads
            && record.clade_percent_classified >= thresholds.min_prop_unmapped_reads
        {
            // If oncogenic_only is true, don't log them even if they pass our thresholds
            if thresholds.oncogenic_only & !is_oncogenic_microbe {
                n_non_oncogenics_excluded += 1;
                continue;
            }

            // Write to our output file
            wtr.serialize(KrakenHit {
                taxid: &record.taxid,
                rank: &record.rank,
                name: &record.name,
                clade_percent_classified: &record.clade_percent_classified,
                clade_nreads_classified: &record.clade_nreads_classified,
                oncogenic: &is_oncogenic_microbe,
            })
            .context("Failed to write KrakenHit")?;

            n_microbial_hits += 1;
            log::info!(
                "Found {} reads from microbe [{}],  ({:4.1}% of all unmapped reads)",
                record.clade_nreads_classified,
                record.name,
                record.clade_percent_classified
            )
        }

        // println!("{:?}", record);
    }

    if thresholds.oncogenic_only && n_non_oncogenics_excluded > 0 {
        log::info!(
            "Skipped reporting {n_non_oncogenics_excluded} microbes despite kraken read support passing thresholds because they were not in our database of oncogenic microbes."
        )
    }

    log::info!("Found {n_microbial_hits} supected microbial hits{oncogenic_only_text}");

    log::info!(
        "Putative kraken hits written to {:#?}",
        &oncogenic_microbe_counts
    );

    Ok(())
}

#[derive(Debug, serde::Deserialize)]
struct KrakenStdRecords {
    classification_status: String,
    sequence_id: String,
    taxid: u64,
    lca_mapping: String,
}

/// Extract reads matching a specific taxid from a bam
fn extract_reads(path_kout: &PathBuf, taxid: u64, path_bam: &PathBuf) -> Result<(), anyhow::Error> {
    // Check if the .kout file exists
    if !path_kout.exists() {
        panic!(
            "Failed to find standard kraken output (.kout) file: {}",
            path_kout.display()
        );
    }

    // Read .kout file using CSV reader
    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .from_path(path_kout)
        .unwrap_or_else(|_| {
            panic!(
                "Failed to parse kraken std output file: {}",
                path_kout.display()
            )
        });

    for result in rdr.deserialize() {
        let record: KrakenStdRecords =
            result.context("Failed to parse record in kraken std output file")?;

        if record.taxid == taxid {}
    }

    Ok(())
}
