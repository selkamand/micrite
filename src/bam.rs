// Take a path to a bam file
use core::str;
use rust_htslib::bam::{self, record::Aux, FetchDefinition, Read};
use rust_htslib::errors::Error;
use std::io::Write;
use std::path::Path;

use crate::kraken::KrakenConfig;

pub fn bam2microbes(bam: &str, outdir: &str, config_kraken: &KrakenConfig) {
    //Filepaths
    let bam_path = std::path::Path::new(bam);
    assert!(
        bam_path.exists(),
        "Could not find BAM file [{}]",
        bam_path.to_str().unwrap()
    );
    let bam_prefix = bam_path
        .file_stem()
        .expect("failed to extract file stem")
        .to_str()
        .expect("Failed to convert bam file stem into prefix");

    let unmapped_fasta = format!("{outdir}/{bam_prefix}.fasta");
    // Create working directory
    std::fs::create_dir_all(outdir).expect("Failed to create output directory");

    // Collect unmapped reads into FASTQAformat
    bam2unmappedreads(bam, unmapped_fasta.as_str(), 50, 17.0);
    eprintln!("Created fasta file of unmapped reads at {unmapped_fasta}");

    // Run Kraken
    let kraken_paths = crate::kraken::run_kraken(unmapped_fasta.clone().into(), config_kraken);

    // Identify Kraken Hits
    crate::kraken::identify_kraken_hits_from_kreport(
        kraken_paths,
        &config_kraken.kraken_hit_thresholds,
    );

    // Delete unmapped fastqs
    eprintln!("Removing unmapped read file");
    std::fs::remove_file(unmapped_fasta).expect("Failed to delete unmapped reads")

    // Extract microbe specific reads for likely hits
    // crate::kraken::extract_reads_from_microbial_hits
}

// Go from bam to unmapped reads
pub fn bam2unmappedreads(bam_path: &str, fasta_output_path: &str, min_len: usize, min_phred: f64) {
    let microbial_contigs = common_microbial_contigs();

    // Create Bam Reader
    let bam_result = bam::IndexedReader::from_path(bam_path);
    let mut bam = match bam_result {
        Ok(value) => value,
        Err(e) => {
            panic!("An error occurred: {:?}", e);
        }
    };

    // Get Bam Header
    let bam_header = bam.header();
    let contigs: Vec<String> = bam_header
        .target_names()
        .iter()
        .map(|t| std::str::from_utf8(t).unwrap().to_string())
        .collect();
    // Braces set to end mutable borrow of bam.header()

    // eprintln!("Bam has the following contigs: {:#?}", contigs);
    let observed_microbial_contigs: Vec<String> = contigs
        .iter()
        .filter(|c| microbial_contigs.contains(c))
        .cloned()
        .collect();

    // Check if we found any microbial contigs
    if !observed_microbial_contigs.is_empty() {
        eprintln!(
            "Found {} contigs in bam that are probably microbial: [{}]",
            observed_microbial_contigs.len(),
            observed_microbial_contigs.join(",")
        )
    }

    // Grab BAM Summary Stats
    let idxstats = bam.index_stats().expect("Failed to get index stats");
    let total_reads: u64 = idxstats.iter().map(|c| c.2 + c.3).sum();
    let total_mapped_reads: u64 = idxstats.iter().map(|c| c.2).sum();
    let total_unmapped_reads: u64 = idxstats.iter().map(|c| c.3).sum();
    eprintln!("BAM-level summary:");
    eprintln!("\ttotal depth (number of reads): [{}]", total_reads);
    eprintln!("\ttotal mapped reads: [{}]", total_mapped_reads);
    eprintln!("\ttotal unmapped reads: [{}]", total_unmapped_reads);
    // Write Bam Summary Stats
    let outdir = Path::new(fasta_output_path)
        .parent()
        .unwrap_or(Path::new("."))
        .to_str()
        .unwrap();

    let stem = Path::new(fasta_output_path)
        .file_stem()
        .unwrap()
        .to_str()
        .unwrap();

    let mut summary_writer = std::fs::File::create(format!("{outdir}/{stem}.bam_summary.txt"))
        .expect("failed to open connection to bam summary stats file");
    writeln!(
        summary_writer,
        "total depth (number of reads)\t{}",
        total_reads
    )
    .expect("Bam summary write failed");
    writeln!(summary_writer, "total mapped reads\t{}", total_mapped_reads)
        .expect("Bam summary write failed");
    writeln!(
        summary_writer,
        "total unmapped reads\t{}",
        total_unmapped_reads
    )
    .expect("Bam summary write failed");

    // Fetch Just the Unmapped reads (based on unmapped flag)
    // Note that some aligners may not set unmapped flag properly
    // (e.g. sometimes if mate read maps the paired unmapped flag is not set).
    // Since the only way to get a complete set of unmapped reads is to manually
    // look through cigar strings of every read, we're going to assume
    // upstream aligners do the right thing.
    bam.fetch(FetchDefinition::Unmapped)
        .expect("Failed to fetch unmapped reads from bam");

    // Open the output FASTA file
    let mut fasta_writer = std::fs::File::create(fasta_output_path)
        .expect("fasta file to output unmapped reads could not be created");

    // Iterate through Unmapped reads and Save to FASTA if they're good quality
    let mut unmapped_good_quality_sequences: u64 = 0;
    let mut unmapped_counter: u64 = 0;
    for r in bam.records() {
        let record = r.unwrap_or_else(|err| panic!("Failed to read bam record: {:?}", err));
        let bam_record = parse_record(&record);
        unmapped_counter += 1;
        // Write to the FASTA file in the correct format
        if is_good_quality_sequence(&bam_record, 50, 17.0, 2) {
            unmapped_good_quality_sequences += 1;
            writeln!(
                fasta_writer,
                ">{}\n{}",
                bam_record.qname, bam_record.sequence
            )
            .expect("Failed to write unmapped read to FASTA file");
        }
    }
    eprintln!("Unmapped Read Summary: ");
    eprintln!("\ttotal unmapped reads: [{}]", unmapped_counter);
    eprintln!(
        "\tgood quality sequences: [{}]",
        unmapped_good_quality_sequences
    );

    // TODO: iterate through any contigs matching known microbial contigs and write mapped reads
    for contig_name in observed_microbial_contigs {
        bam.fetch(&contig_name)
            .expect("Error fetching bam sequences from specific contigs");

        let mut nreads: u64 = 0;
        let mut nreads_mapped: u64 = 0;
        let mut nreads_good_sequence: u64 = 0;
        let mut nreads_good_alignment: u64 = 0;
        for r in bam.records() {
            let record = r.unwrap_or_else(|err| panic!("Failed to read bam record: {:?}", err));
            let bam_record = parse_record(&record);

            nreads += 1;

            if (!record.is_unmapped()) {
                nreads_mapped += 1
            }

            // Write good quality sequences mapped to microbial contigs to the fasta file
            if !record.is_unmapped() & is_good_quality_sequence(&bam_record, 50, 17.0, 2) {
                nreads_good_sequence += 1;
                writeln!(
                    fasta_writer,
                    ">{}\n{}",
                    bam_record.qname, bam_record.sequence
                )
                .expect("Failed to write unmapped read to FASTA file");
            }

            // Count Number of Good Quality Alignments
            // TODO: MAke alignment scores (AS) sequence length independent (might end up making micrite even more aligner specific though)
            if is_good_quality_alignment(&bam_record, 50, 17.0, 2, 10, 130) {
                nreads_good_alignment += 1
            }
        }
        eprintln!("Microbial Contig Stats: {}", contig_name);
        eprintln!("\ttotal reads mapped: [{}]", nreads_mapped);
        eprintln!(
            "\tgood quality alignments mapped: [{}]",
            nreads_good_alignment
        );
        eprintln!(
            "\tgood quality sequences mapped: [{}]",
            nreads_good_sequence
        );
        writeln!(
            summary_writer,
            "Contig [{}] good quality alignments\t{}",
            contig_name, nreads_good_alignment
        )
        .expect("Failed write");
    }
}

// A custom struct that adds a couple of key properties to bam::record
struct BamRecordEnriched<'a> {
    record: &'a rust_htslib::bam::Record,
    qname: &'a str,
    sequence: String,
    alignment_score: i32,
}

fn get_as_tag(record: &bam::Record) -> Option<i32> {
    match record.aux(b"AS") {
        Ok(Aux::I8(value)) => Some(value as i32),
        Ok(Aux::U8(value)) => Some(value as i32),
        Ok(Aux::I16(value)) => Some(value as i32),
        Ok(Aux::U16(value)) => Some(value as i32),
        Ok(Aux::I32(value)) => Some(value),
        Ok(Aux::U32(value)) => Some(value as i32),
        Ok(_) => None, // The AS tag exists but is of an unexpected type
        Err(Error::BamAuxTagNotFound) => None, // AS tag not found
        Err(e) => {
            // Handle other potential errors
            eprintln!("Error retrieving AS tag: {}", e);
            None
        }
    }
}

fn parse_record(record: &bam::Record) -> BamRecordEnriched {
    // Run computationally intensive checks
    let seq = record.seq().as_bytes();
    let sequence: String = seq.iter().map(|&b| b as char).collect();
    let qname = str::from_utf8(record.qname()).expect("Failed to parse qname to string slice");
    let alignment_score = get_as_tag(record).unwrap_or(0);

    BamRecordEnriched {
        record,
        qname,
        sequence,
        alignment_score,
    }
}

/// Check whether a bam sequence is considered 'good quality'.
///
/// A good quality *sequence* is likely to be a real biological
/// sequence that should be fed into kraken downstream for read classification.
/// Note a good quality sequence is not necessarily a good quality 'alignment'
///
/// A good quality sequence has the following properties
/// 1. Reasonable length (>`min_len``)
/// 2. Good Average Phred Scores (>=`min_phred`)
/// 3. Contains very few ambiguous/masked nucleotides (Number of Ns < `max_n`)
/// 4. Is not a PCR duplicate or flagged as 'is_quality_check_failed'
/// 5. Has a reasonable sequence complexity (No homopolymer reads) (not yet implemented)
///
fn is_good_quality_sequence(
    record: &BamRecordEnriched,
    min_len: usize,
    min_phred: f64,
    max_n: usize,
) -> bool {
    // Start with the quick checks

    if record.record.is_quality_check_failed()
        | record.record.is_duplicate()
        | (record.record.seq_len() < min_len)
    {
        return false;
    }

    // Run computationally intensive checks
    // Ambiguous bases (N)
    let has_ambiguous_bases: bool = seq_ambiguous(&record.sequence, max_n);

    // Average Quality
    let qual = record.record.qual();
    let qual_average = calculate_average_phred(qual);

    if has_ambiguous_bases | (qual_average < min_phred) {
        return false;
    }

    // TODO: Add a check based on sequence complexity

    return true;
}

/// Is the alignment convincing
fn is_good_quality_alignment(
    record: &BamRecordEnriched,
    min_len: usize,
    min_phred: f64,
    max_n: usize,
    min_mapq: u8,
    min_alignment_score: i32,
) -> bool {
    // CHeck if sequence is good quality
    let good_qual_sequence = is_good_quality_sequence(record, min_len, min_phred, max_n);
    if !good_qual_sequence {
        return false;
    }

    // Check if Alignment is good quality
    //TODO: add an aditional check on absolute mapping quality between seq and ref (Maybe using AS tag)
    !record.record.is_secondary()
        & !record.record.is_quality_check_failed()
        & !record.record.is_unmapped()
        & (record.record.mapq() > min_mapq)
        // Alignment Score 
        & (record.alignment_score > min_alignment_score)
}

/// Check how many Ns in a string, and if greater than 'maxNs' return FALSE
fn seq_ambiguous(seq: &str, max_n: usize) -> bool {
    let number_of_ns = seq.chars().filter(|c| *c == 'N').count();
    number_of_ns > max_n
}

fn calculate_average_phred(qual_scores: &[u8]) -> f64 {
    let total: u32 = qual_scores.iter().map(|&score| score as u32).sum();
    let count = qual_scores.len();

    if count > 0 {
        total as f64 / count as f64
    } else {
        0.0
    }
}

struct SeqClassification {
    ambiguous: bool,
    low_complexity: bool,
}
#[derive(Debug, serde::Deserialize)]
struct MicrobialContigRecords {
    taxid: String,
    common_name: String,
    contigs: String,
}
struct Contig {
    contig: String,
    taxid: String,
    species: String,
}

/// A collection of microbial contigs.
/// Use the `contains` method to see if a particular contig name is in the list
pub struct MicrobialContigs {
    contigs: Vec<Contig>,
}

impl MicrobialContigs {
    // Check if InterestingContigs contain a particular contig name
    fn contains(&self, contig_name: &str) -> bool {
        let contigs_in_set: Vec<&str> = self.contigs.iter().map(|c| c.contig.as_str()).collect();
        contigs_in_set.contains(&contig_name)
    }

    // If Taxid
    fn contig_to_species(&self, contig_name: &str) -> Option<&str> {
        let species = self
            .contigs
            .iter()
            .filter(|c| c.contig.as_str() == contig_name)
            .map(|c| c.species.as_str())
            .next();

        species
    }
}

pub fn common_microbial_contigs() -> MicrobialContigs {
    MicrobialContigs {
        contigs: vec![
            //EBV
            Contig {
                contig: "chrEBV".to_string(),
                taxid: "10376".to_string(),
                species: "EBV".to_string(),
            },
            Contig {
                contig: "NC_009334".to_string(),
                taxid: "10376".to_string(),
                species: "EBV".to_string(),
            },
            Contig {
                contig: "NC_007605".to_string(),
                taxid: "10376".to_string(),
                species: "EBV".to_string(),
            },
            //HHV6B
            Contig {
                contig: "NC_000898".to_string(),
                taxid: "10376".to_string(),
                species: "HHV6B".to_string(),
            },
        ],
    }
}

#[cfg(test)]
mod tests {

    #[test]
    fn microbial_contigs() {
        let microcontigs = crate::bam::common_microbial_contigs();
        assert!(microcontigs.contains("NC_007605"));
        assert_eq!(microcontigs.contig_to_species("NC_007605").unwrap(), "EBV");
        assert_eq!(
            microcontigs.contig_to_species("NC_000898").unwrap(),
            "HHV6B"
        );
        assert!(microcontigs.contig_to_species("ADAWD").is_none());
    }
}
