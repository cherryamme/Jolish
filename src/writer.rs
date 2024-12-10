// use std::collections::HashMap;
use flate2::write::GzEncoder;
use flate2::Compression;
use log::debug;
use log::info;
use noodles::bam::record;
use std::fs::File;
use std::hash::Hash;
use std::io::Write;
// use std::io::Result;
use std::path::Path;
// use bio::io::fastq::Writer;
use crate::bam::ReadPart;
use crate::bam::ReadPartChunk;
use crate::bam::ReadpartsNumDictValue;
use crate::fastq::ReadChunk;
use crate::fastq::ReadInfo;
use std::fs::create_dir_all;
use std::io::BufWriter;
// use std::thread;
use flume::{unbounded, Receiver, Sender};
use std::collections::HashMap;
use std::thread;
use std::time::Instant;

use std::sync::{Arc, Mutex};

// 创建一个函数，输入sequence序列，对碱基序列进行反向互补，返回一个反向互补的字符串
fn reverse_complement(seq: &str) -> String {
    seq.chars().rev().map(complement).collect()
}

fn complement(base: char) -> char {
    match base {
        'A' => 'T',
        'T' => 'A',
        'C' => 'G',
        'G' => 'C',
        _ => base,
    }
}


// 创建一个函数，输入read_parts,orient,将read_parts按照read_rank排序，返回一个sequence和quality字符串
// pub fn assemble_read_parts(read_parts: &[&ReadPart], orient: &str) -> (String, String) {
//     let mut record_seq = String::new();
//     let mut record_qual = String::new();

//     for rp in read_parts {
//         let (seq, qual) = if rp.correct {
//             (&rp.correct_seq, &rp.correct_qual)
//         } else {
//             (&rp.read_seq, &rp.read_qual)
//         };
//         debug!("seq is {:?}, qual is {:?}", seq, qual);
//         record_seq.push_str(std::str::from_utf8(seq).expect("Not a valid UTF-8 sequence"));
//         record_qual.push_str(std::str::from_utf8(qual).expect("Not a valid UTF-8 sequence"));
//     }

//     let seq = if orient == "+" {
//         record_seq
//     } else {
//         reverse_complement(&record_seq)
//     };

//     let qual = if orient == "+" {
//         record_qual
//     } else {
//         record_qual.chars().rev().collect()
//     };
//     debug!("seq is {:?}, qual is {:?}", seq, qual);
//     (seq, qual)
// }



// 创建一个函数，输入record_info，read_parts，将read_parts按照read_rank排序，返回一个record_info
// pub fn assemble_read(record_info: String, read_parts: &Vec<ReadPart>) -> String {

// }
pub fn assemble_read_parts(read_parts: &[&ReadPart], orient: &str) -> (String, String) {
    let mut record_seq = String::new();
    let mut record_qual = Vec::new();

    for rp in read_parts {
        let (seq, qual) = if rp.correct {
            (&rp.correct_seq, &rp.correct_qual)
        } else {
            (&rp.read_seq, &rp.read_qual)
        };
        // Remove the first and last eight bases
        let trimmed_seq = &seq[8..seq.len() - 8];
        let trimmed_qual = &qual[8..qual.len() - 8];

        record_seq.push_str(std::str::from_utf8(trimmed_seq).expect("Not a valid UTF-8 sequence"));
        // Add 33 to each quality value since they've been subtracted by 33
        record_qual.extend(trimmed_qual.iter().map(|&q| q + 33));
    
    }

    let seq = if orient == "+" {
        record_seq
    } else {
        reverse_complement(&record_seq)
    };

    let qual = if orient == "+" {
        String::from_utf8(record_qual).expect("Not a valid UTF-8 sequence")
    } else {
        let reversed_qual: Vec<u8> = record_qual.into_iter().rev().collect();
        String::from_utf8(reversed_qual).expect("Not a valid UTF-8 sequence")
    };
    // debug!("seq is {:?}, qual is {:?}", seq, qual);
    (seq, qual)
}




pub fn writer_receiver_bam(
    rrx: Receiver<ReadPartChunk>,
    readparts_num_dict: &Mutex<HashMap<String, ReadpartsNumDictValue>>,
    outfile: &str,
) {
    // 创建一个HashMap，key为read_id，value为readpart
    let mut writer_fastq_dict = HashMap::new();
    let start_time = Instant::now();
    let filepath = Path::new(outfile);
    let filedir = &filepath.parent().unwrap();
    create_dir_all(&filedir).expect("fail to create output directory");
    let file = File::create(&filepath).expect("fail to create output fq.gz");
    let encoder = GzEncoder::new(file, Compression::default());
    // 创建一个BufWriter，用于写入文件
    let mut writer = BufWriter::with_capacity(1_000_000, encoder);
    let mut readpart_count = 0;
    for readchunk in rrx {
        for readpart in readchunk {
            readpart_count += 1;
            let read_id = readpart.read_id.clone();
            // info!("loading readpart {:?}",readpart);
            // 将readpart按照read_id存入writer_fastq_dict
            writer_fastq_dict
                .entry(read_id.clone())
                .or_insert(Vec::new())
                .push(readpart.clone());

            let readparts_num_map = readparts_num_dict
                .lock()
                .expect("Failed to lock readparts_num_dict");
            
            if let Some(readpartsnumdictvalue) = readparts_num_map.get(&read_id) {
                let read_parts = writer_fastq_dict.get(&read_id).unwrap();
                // debug!(
                //     "read_parts_num is {:?} read_parts: {:?}",
                //     read_parts_num,
                //     read_parts.len()
                // );
                // 如果writer_fastq_dict中这个readid的reads数量达到了readparts_num_dict中的chunk_num的数量，并且ReadpartsNumDictValue中的record_num==record_count，则写入文件
                if read_parts.len() == readpartsnumdictvalue.chunk_num
                    && readpartsnumdictvalue.record_num == readpartsnumdictvalue.record_count
                {
                    // Sort read parts by read_rank
                    // TODO reads 质量值 待修复
                    // TODO 创建一个函数，输入record_info，read_parts，将read_parts按照record_info分组处理，返回一个record_str

                    let mut sorted_parts: Vec<_> = read_parts.iter().collect();
                    sorted_parts.sort_by_key(|rp| rp.read_rank);
                    let (record_seq,record_qual) = assemble_read_parts(&sorted_parts, &readpartsnumdictvalue.orient);
                    let record_str =
                        format!("@{}\n{}\n+\n{}\n", read_id, record_seq, record_qual);
                    // info!("record_str is {:?}", record_str);
                    write!(writer, "{}", record_str).expect("Failed to write to writer");
                    // Optionally clear the entry if no longer needed
                    writer_fastq_dict.remove(&read_id);
                }
            }
        }
    }
    let elapsed_time = start_time.elapsed();
    info!(
        "{} readpart writer to fastq {}.\n Time elapsed: {:.4?}",
        readpart_count, outfile, elapsed_time
    )
}



pub fn writer_receiver(rrx: Receiver<ReadChunk>, outfile: &str) {
    let start_time = Instant::now();
    let filepath = Path::new(outfile);
    let filedir = &filepath.parent().unwrap();
    create_dir_all(&filedir).expect("fail to create output directory");
    let file = File::create(&filepath).expect("fail to create output fq.gz");
    let encoder = GzEncoder::new(file, Compression::default());
    let mut writer = BufWriter::with_capacity(1_000_000, encoder);
    for readchunk in rrx {
        for readinfo in readchunk.reads {
            let id = readinfo.record.id();
            let seq =
                std::str::from_utf8(&readinfo.correct_seq).expect("Not a valid UTF-8 sequence");
            let qual =
                std::str::from_utf8(&readinfo.correct_qual).expect("Not a valid UTF-8 sequence");
            let record_str = format!("@{}\n{}\n+\n{}\n", id, seq, qual);
            write!(writer, "{}", record_str).unwrap();
        }
    }
    let elapsed_time = start_time.elapsed();
    info!(
        "writer to fastq {}.\n Time elapsed: {:.4?}",
        outfile, elapsed_time
    )
}
