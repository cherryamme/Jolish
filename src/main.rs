#![allow(warnings)]
use log::{info,debug};
use clap::Parser;
mod args;
pub mod poa_test;
pub mod poa;
pub mod correcter;
pub mod fastq;
pub mod writer;
pub mod fastq_handler;
pub mod bam_handler;
pub mod bam;


fn main() {
    std::env::set_var("RUST_BACKTRACE", "1");
    std::env::set_var("RUST_LOG", "debug");
    pretty_env_logger::init();
    let comands: Vec<String> = std::env::args().collect();
    let mut args = args::Args::parse();
    info!("Run Command: {:?}", comands);
    info!("Args: {:?}", args);
    // // 判断inputs是bam还是fastq，使用不同的handler
    // args.inputs.push("/home/jiangchen/project/Jolish/benchmark2/benchmark_1000.fq".to_string());
    // args.outfile = "/home/jiangchen/project/Jolish/benchmark2/benchmark_1000_correct.fq.gz".to_string();
    if args.inputs[0].ends_with(".bam") {
        let (rrx, readparts_num_dict) = bam::spawn_bam_reader(args.inputs, args.chunk_size, args.chunk_window, args.limit);
        let srx = bam_handler::correcter_receiver(rrx, args.threads, args.correct_ratio, args.correct_fix_ratio);
        // for readchunk in srx.iter(){
        //     for read in readchunk.iter(){
        //         info!("read id is {:?}",read);
        //     }
        // }
        // sleep 20s
        std::thread::sleep(std::time::Duration::from_secs(1));
        // debug!("{:?}",readparts_num_dict);
        writer::writer_receiver_bam(srx, &readparts_num_dict, &args.outfile);
        return;
    }else if args.inputs[0].ends_with(".fq") || args.inputs[0].ends_with(".fastq") || args.inputs[0].ends_with(".gz") {
        let rrx: flume::Receiver<fastq::ReadChunk> = fastq::spawn_reader(args.inputs, args.chunk_size, args.split_param, args.type_index, args.orient_index, args.limit);
        let srx = fastq_handler::correcter_receiver(rrx, args.threads, args.correct_ratio, args.correct_fix_ratio);
        writer::writer_receiver(srx, &args.outfile);
        return;
    }
}


#[test]
fn test_rrx(){
    std::env::set_var("RUST_LOG", "debug");
    pretty_env_logger::init();
    let rrx: flume::Receiver<fastq::ReadChunk> = fastq::spawn_reader(vec!["/home/jiangchen/project/Jolish/benchmark/test.fq.gz".to_string()], 10, '%', 1, 3, Some(200));
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
