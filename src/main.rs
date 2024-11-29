use log::{info,debug};
use clap::Parser;
mod args;
pub mod poa_test;
pub mod poa;
pub mod correcter;
pub mod fastq;
pub mod writer;
pub mod correct_handler;



fn main() {
    std::env::set_var("RUST_LOG", "info");
    pretty_env_logger::init();
    let comands: Vec<String> = std::env::args().collect();
    let args = args::Args::parse();
    info!("Run Command: {:?}", comands);
    info!("Args: {:?}", args);
    let rrx: flume::Receiver<fastq::ReadChunk> = fastq::spawn_reader(args.inputs, args.chunk_size, args.split_param, args.type_index, args.orient_index, args.max_reads);
    // let rrx: flume::Receiver<fastq::ReadChunk> = fastq::spawn_reader(vec!["/home/jiangchen/project/Jolish/benchmark/benchmark_500.fq".to_string()], args.chunk_size, args.split_param, args.type_index, args.orient_index, args.max_reads);
    let srx = correct_handler::correcter_receiver(rrx, args.threads, args.chunk_size, args.correct_ratio, args.correct_fix_ratio);
    writer::writer_receiver(srx, &args.outfile);
}


#[test]
fn test_rrx(){
    std::env::set_var("RUST_LOG", "debug");
    pretty_env_logger::init();
    let rrx: flume::Receiver<fastq::ReadChunk> = fastq::spawn_reader(vec!["/home/jiangchen/project/Jolish/test.fq".to_string()], 10, '%', 1, 3, Some(200));
    for mut chunk in rrx.iter() {
        debug!("chunk length is {}",chunk.reads.len());
        debug!("key is {:?}",chunk.info());
        for read in chunk.iter_mut() {
            debug!("read id is {}",read.info());
            debug!("read seq is {:?}",read.record.seq());
        }
    }
    
}

// #[test]
// fn test_correct_handler(){
//     std::env::set_var("RUST_LOG", "debug");
//     pretty_env_logger::init();
//     let rrx: flume::Receiver<fastq::ReadChunk> = fastq::spawn_reader(vec!["/home/jiangchen/project/Jolish/test.fq".to_string()], 10, '%', 1, 3);
//     for readchunk in rrx.iter(){
//         let mut ch =  correct_handler::CorrectHandler::new(readchunk);
//         ch.info();
//         // ch.correct_reads();
//     }
    
// }
