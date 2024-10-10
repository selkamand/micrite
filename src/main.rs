fn main() {
    // Screen BAM for microbial reads using a kraken2 database
    // micrite::bam2unmappedreads(bam_path, bam_output_path);
    // bam = "inst/"

    let fasta = micrite::bam::bam2microbes(
        "testfiles/humanGRCh38_9000_ebv_1000_hpv16_1000_hpylori_1000.grch38_noalt.bam",
        "outdir",
    );

    let config = micrite::kraken::KrakenConfig {
        krakendb: std::path::PathBuf::from("~/databases/kraken2/k2_standard_08gb_20240605"),
        threads: 8,
        confidence: "0.01".to_string(),
        outdir: "outdir".to_string(),
    };
    micrite::kraken::run_kraken(fasta, config);
}
