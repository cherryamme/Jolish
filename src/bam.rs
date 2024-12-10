use flate2::read;
use flume::{unbounded, Receiver, Sender};
use log::{debug, error, info};
use noodles::bam::{self, record};
use noodles::sam::alignment::record::cigar::op::Kind;
use noodles::sam::alignment::record::data::field::Tag;
use noodles::sam::alignment::record::data::field::Value;
use petgraph::graph::NodeIndex;
use std::collections::HashMap;
use std::{clone, io};
use std::sync::{Arc, Mutex};
use std::time::Instant;
use std::{
    fs::File,
    io::{BufReader, Read},
    path::PathBuf,
};

pub type ReadPartChunk = Vec<ReadPart>;
const SA: [u8; 2] = [b'S', b'A'];
#[derive(Debug)]
pub struct ReadpartsNumDictValue {
    pub orient: String,
    pub record_info: Vec<RecordInfo>,
    pub record_num: usize,
    pub record_count: usize,
    pub chunk_num: usize,
}

#[derive(Debug,Clone)]
pub struct RecordInfo {
    pub record_id: String,
    pub record_order: usize,
    pub orient: String,
}

#[derive(Debug)]
pub struct ReadParts_idx {
    pub chrom: usize,
    pub windows_num: usize,
    pub read_order: usize,
    pub start: usize,
    pub end: usize,
    pub ref_pos: usize,
}
#[derive(Debug, Clone)]
pub struct ReadPart {
    pub read_id: String,
    pub record_id: String,
    pub read_rank: usize,
    pub read_seq: Vec<u8>,
    pub read_qual: Vec<u8>,
    pub seq_node_index: Vec<NodeIndex<usize>>,
    pub correct_seq: Vec<u8>,
    pub correct_qual: Vec<u8>,
    pub correct: bool,
}
impl ReadPart {
    pub fn new(
        read_id: String,
        record_id: String,
        read_seq: &Vec<u8>,
        read_qual: &Vec<u8>,
        readparts_idx: &ReadParts_idx,
    ) -> ReadPart {
        let correct = readparts_idx.end - readparts_idx.start >= 70;
        let mut seq = read_seq[readparts_idx.start..readparts_idx.end].to_vec();
        let mut qual = read_qual[readparts_idx.start..readparts_idx.end].to_vec();
        // Add "ATGCATGC" to the beginning and end of the sequence
        let prefix_suffix = b"ATGCATGC".to_vec();
        seq.splice(0..0, prefix_suffix.clone());
        seq.extend(prefix_suffix.clone());

        // Add Q30 quality values to the beginning and end of the quality vector
        let q30 = vec![30; 8];
        qual.splice(0..0, q30.clone());
        qual.extend(q30.clone());

        ReadPart {
            read_id: read_id,
            record_id: record_id,
            read_rank: readparts_idx.read_order,
            read_seq: seq,
            read_qual: qual,
            seq_node_index: Vec::new(),
            correct_seq: Vec::new(),
            correct_qual: Vec::new(),
            correct: correct,
        }
    }
    pub fn need_correct(&self) -> bool {
        self.correct
    }
}

fn get_readparts_idx(record: &bam::Record, chunk_window: usize) -> Vec<(ReadParts_idx)> {
    let cigar = record.cigar();
    let reference_id = record.reference_sequence_id().unwrap().unwrap();
    let mut readparts_idx = Vec::new();
    // debug!("{:?}",cigar);
    let start = record.alignment_start().unwrap().unwrap().get();
    let mut windows_num = start / chunk_window;
    let windows_left = start % chunk_window;
    let mut read_index = 0;
    let mut need_bases = chunk_window - windows_left;
    let mut start_pos = 0;
    let mut read_order = 0;
    let mut ref_pos = start - windows_left;
    // debug!(
    //     "start:{} windows_num:{} windows_left:{}",
    //     start, windows_num, windows_left
    // );
    for c in cigar.iter() {
        let op = c.unwrap();
        // debug!("kind:{:?} len:{:?}",op.kind(),op.len());
        match op.kind() {
            Kind::Match => {
                let mut op = op.len();
                while op >= need_bases {
                    // info!(
                    //     "Match:windows_num:{} windows_left:{} read_index:{} need_bases:{}",
                    //     windows_num, windows_left, read_index, need_bases
                    // );
                    windows_num += 1;
                    op -= need_bases;
                    read_index += need_bases;
                    need_bases = chunk_window;
                    // send readparts
                    // start_pos-read_index
                    // debug!("Match: Send readparts start:{} end:{}",start_pos,read_index);
                    let readpart = ReadParts_idx {
                        chrom: reference_id,
                        windows_num: windows_num,
                        read_order: read_order,
                        start: start_pos,
                        end: read_index,
                        ref_pos: ref_pos,
                    };
                    readparts_idx.push(readpart);
                    start_pos = read_index;
                    read_order += 1;
                    ref_pos += chunk_window; // Increment ref_pos by 100
                }
                need_bases -= op;
                read_index += op;
            }
            Kind::Insertion => {
                read_index += op.len();
            }
            Kind::Deletion => {
                let mut op = op.len();
                while op >= need_bases {
                    // debug!("Deletion:windows_num:{} windows_left:{} read_index:{} need_bases:{}",windows_num,windows_left,read_index,need_bases);
                    windows_num += 1;
                    op -= need_bases;
                    need_bases = chunk_window;
                    // send readparts
                    // start_pos-read_index
                    if start_pos != read_index {
                        // debug!("Deletion: Send readparts start:{} end:{}",start_pos,read_index);
                        let readpart = ReadParts_idx {
                            chrom: reference_id,
                            windows_num: windows_num,
                            read_order: read_order,
                            start: start_pos,
                            end: read_index,
                            ref_pos: ref_pos,
                        };
                        readparts_idx.push(readpart);
                        start_pos = read_index;
                        read_order += 1;
                        ref_pos += chunk_window;
                    }
                }
            }
            Kind::SoftClip => {
                // debug!("SoftClip:windows_num:{} windows_left:{} read_index:{} need_bases:{}",windows_num,windows_left,read_index,need_bases);
                let op = op.len();
                if start_pos == 0 {
                    start_pos += op;
                    read_index += op;
                    read_order += (op / 200) + 1;
                } else {
                    // send readparts
                    // start_pos-read_index
                    // debug!("SoftClip: Send readparts start:{} end:{}",start_pos,read_index);
                    let readpart = ReadParts_idx {
                        chrom: reference_id,
                        windows_num: windows_num,
                        read_order: read_order,
                        start: start_pos,
                        end: read_index,
                        ref_pos: ref_pos,
                    };
                    readparts_idx.push(readpart);
                    read_order += 1;
                }
            }
            Kind::HardClip => {
                error!("HardClip");
            }
            Kind::Pad => {
                error!("Pad");
            }
            Kind::Skip => {
                error!("Skip");
            }
            Kind::SequenceMatch => {
                error!("SequenceMatch");
            }
            Kind::SequenceMismatch => {
                error!("SequenceMismatch");
            }
        }
    }
    // 结尾如果还有剩余，需要发送
    if start_pos != read_index {
        // debug!("End: Send readparts start:{} end:{}",start_pos,read_index);
        let readpart = ReadParts_idx {
            chrom: reference_id,
            windows_num: windows_num,
            read_order: read_order,
            start: start_pos,
            end: read_index,
            ref_pos: ref_pos,
        };
        readparts_idx.push(readpart);
    }
    readparts_idx
}

fn parse_satag_to_count(value: Option<Result<Value<'_>, io::Error>>) -> usize {
    match value {
        None => 1,
        Some(Err(e)) => 1, // Convert the error to a string
        Some(Ok(v)) => {
            let content = format!("{:?}", v);
            content.split(";").count()
        } // Use the Debug trait to convert the value to a string
    }
}

pub fn spawn_bam_reader(
    files: Vec<String>,
    chunk_size: usize,
    chunk_window: usize,
    max_reads: Option<usize>,
) -> (
    Receiver<ReadPartChunk>,
    Arc<Mutex<HashMap<String, ReadpartsNumDictValue>>>,
) {
    let (rtx, rrx) = unbounded();
    let readparts_num_dict = Arc::new(Mutex::new(HashMap::new()));
    // Create a mutable reference to `readparts_num_dict`
    let readparts_num_dict_clone = Arc::clone(&readparts_num_dict); // Clone the Arc for use in the thread
    std::thread::spawn(move || {
        let start_time = Instant::now();
        let mut readparts_dict = HashMap::new();
        let mut record_count = 0;
        let mut readparts_count = 0;
        for file in files {
            // 加载bam文件
            let mut reader = bam::io::reader::Builder::default()
                .build_from_path(&file)
                .unwrap();
            let header = reader.read_header().unwrap();
            // log日志，打印bam文件的header信息
            info!("Process bam: {:?}\nHeader: {:?}", file, header);
            use noodles::sam::alignment::record::data::field::Value;
            // 循环处理bam文件中的每一行
            for result in reader.records() {
                let record = result.unwrap();
                // 如果record是unmapped的，跳过
                if record.flags().is_unmapped() {
                    continue;
                }
                // 处理的record数量加1
                record_count += 1;
                // 读取record是正向还是反向
                let orient = if record.flags().is_reverse_complemented() {
                    "reverse"
                } else {
                    "forward"
                };
                // 解析record的信息
                let read_id: String = record.name().unwrap().to_string();
                let read_seq: Vec<u8> = record.sequence().iter().collect();
                let read_qual: Vec<u8> = record.quality_scores().as_ref().to_vec();
                // 计算record的分窗信息
                let readparts_idx = get_readparts_idx(&record, chunk_window);
                if read_id == "TB2000CC3C-202406061734310_20240606183107_00177_000021521_11.46%4.2-F_3.7-R%00001.split.fastq.gz%-" {
                    debug!("readparts_idx : {:?}", readparts_idx);
                }
                // 计算record的supplementary_alignment的数量，该reads的比对数量
                let record_data = record.data();
                let sa_tag: Option<Result<Value<'_>, std::io::Error>> = record_data.get(&SA);
                let sa_tag_count = parse_satag_to_count(sa_tag);
                // 如果比对上多个地方，直接跳过，进行下一个record的校正
                if sa_tag_count > 1 {
                    continue;
                }
                // TODO 反向互补待增加
                // 解析当前record的信息，生成record_id，并生成record_dict，内部保留record_order和orient信息
                let record_order = readparts_idx[0].read_order;
                let record_id = format!("{}_{}", read_id, record_order);

                // 生成record_info结构体，内部包含record_id, record_order, orient
                let record_info = RecordInfo {
                    record_id: record_id.clone(),
                    record_order: record_order,
                    orient: orient.to_string(),
                };
                // debug!(
                //     "record_id: {:?} readparts_idx: {:?} sa_tag_count: {:?}",
                //     record_id, readparts_idx, sa_tag_count
                // );
                // lock然后更新readparts_num_dict
                // 给readparts_num_dict内的readid的值增加 record_info, record_count, chunk_num
                let mut dict = readparts_num_dict_clone.lock().unwrap();
                dict.entry(read_id.clone())
                    .and_modify(|value: &mut ReadpartsNumDictValue| {
                        // 根据需要更新其他字段
                        value.record_info.push(record_info.clone());
                        value.record_count += 1; // 假设每次调用增加一次记录计数
                        value.chunk_num += readparts_idx.len(); // 假设每次调用增加一个块
                    })
                    .or_insert(ReadpartsNumDictValue {
                        record_info: vec![record_info],
                        orient: orient.to_string(),
                        record_num: sa_tag_count,
                        record_count: 1,                // 初始化为1
                        chunk_num: readparts_idx.len(), // 初始化为1
                    });
                // 根据record的分窗信息，生成readpart，并将readpart发送到rtx
                for readpart_idx in readparts_idx.iter() {
                    // debug!("readpart_idx: {:?}",readpart_idx);
                    // 分窗内创建readpart结构
                    let readpart =
                        ReadPart::new(read_id.clone(),record_id.clone(), &read_seq, &read_qual, readpart_idx);
                    if read_id == "TB2000CC3C-202406061734310_20240606183107_00177_000021521_11.46%4.2-F_3.7-R%00001.split.fastq.gz%-" {
                        debug!("readpart : {:?}", readpart);
                    }
                    // 将readpart根据染色体号和坐标信息为键合并进字典
                    let dict_key = (readpart_idx.chrom, readpart_idx.windows_num);
                    readparts_dict
                        .entry(dict_key)
                        .or_insert_with(Vec::new)
                        .push(readpart.clone());

                    // if dict_key.0 == 0 && dict_key.1 == 178 {
                    //     debug!("readparts_idx: {:?}",readparts_idx);
                    //     debug!("read_id:{:?} read_seq: {:?}",read_id,std::str::from_utf8(&readpart.read_seq));
                    // }
                    let dict_value = readparts_dict.get_mut(&dict_key).unwrap();
                    // 计算当前字典的键值中需要校正的readpart数量
                    let correct_count = dict_value.iter().filter(|rp| rp.need_correct()).count();
                    // 如果需要校正的readpart数量达到100，发送到rtx进入下游correct channel
                    if correct_count == chunk_size {
                        // debug!(
                        //     "Send dict_value Chunk reads {:?} to correct chanel: {:?}",
                        //     correct_count, dict_key
                        // );
                        let to_send = readparts_dict.remove(&dict_key).unwrap(); // Move the Vec
                        readparts_count += to_send.len();
                        rtx.send(to_send).expect("Error sending");
                    }
                    // info!("readparts_idx: {:?}",readparts_idx);
                    // info!("readparts_dict_len: {:?}",readparts_dict.len());
                    // dbg!(record);
                }
            }
            // 将所有readparts_dict内残留的readpart 发送到correct channel
            for (key, value) in readparts_dict.drain() {
                if !value.is_empty() {
                    // info!(
                    //     "Send dict_value Chunk reads {:?} to correct chanel: {:?}",
                    //     value.len(),
                    //     key
                    // );
                    readparts_count += value.len();
                    rtx.send(value).expect("Error sending");
                }
            }
            let elapsed_time = start_time.elapsed();
            info!("Loading Reads data done! Process record_count: {:?} readpart_count: {:?} Time elapsed: {:.4?}", record_count, readparts_count, elapsed_time);
        }
    });
    (rrx, readparts_num_dict)
}

#[test]
pub fn test_spawn_bam_reader() {
    pretty_env_logger::init();
    let files = vec!["/home/jiangchen/project/Jolish/example/bam_example.fq.bam".to_string()];
    let rrx = spawn_bam_reader(files, 100, 100, Some(200));
    // for readparts in rrx.iter() {
    //     info!("Readparts: {:?}",readparts);
    // }
}

// #[test]
// pub fn test_bam() {
//     let mut readparts_dict = HashMap::new();
//     pretty_env_logger::init();
//     let mut reader = bam::io::reader::Builder::default()
//         .build_from_path("/home/jiangchen/project/Jolish/example/bam_example.fq.bam")
//         .unwrap();
//     let header = reader.read_header().unwrap();
//     info!("{:?}", header);
//     let mut record_count = 0;
//     // 每一行record生成一个fastq中的reads，而不是每一个read生成一个record
//     for result in reader.records() {
//         // debug!("loading_Result: {:?}",result);
//         // let record = ReadInfo::new(result.unwrap(), '%', 1, 3);
//         // debug!("{:?}",record);
//         let record = result.unwrap();
//         if record.flags().is_unmapped() {
//             continue;
//         }
//         if record.flags().is_reverse_complemented() {
//             debug!("reverse_complemented reads");
//         }
//         let read_id: String = record.name().unwrap().to_string();
//         let read_seq: Vec<u8> = record.sequence().iter().collect();
//         let read_qual: Vec<u8> = record.quality_scores().as_ref().to_vec();
//         //
//         let readparts_idx = get_readparts_idx(&record);
//         // debug!("readparts_idx: {:?}",readparts_idx);
//         // 这个record是一个read的所有的需要拆分的信息。
//         // let mut readparts = Vec::new();
//         record_count += 1;
//         for readpart_idx in readparts_idx.iter() {
//             // debug!("readpart_idx: {:?}",readpart_idx);
//             let readpart = ReadPart::new(read_id.clone(),record, &read_seq, &read_qual, readpart_idx);
//             // readparts.push(readpart);
//             let dict_key = (readpart_idx.chrom, readpart_idx.windows_num);
//             if !readparts_dict.contains_key(&dict_key) {
//                 readparts_dict.insert(dict_key, vec![readpart]);
//                 debug!(
//                     "insert new dict_value: {:?}, dict_key: {:?}",
//                     record_count, dict_key
//                 );
//             } else {
//                 let dict_value: &mut Vec<ReadPart> = readparts_dict.get_mut(&dict_key).unwrap();
//                 dict_value.push(readpart);
//                 let correct_count = dict_value.iter().filter(|rp| rp.need_correct()).count();
//                 // debug!("dict_value {:?}: {:?}",record_count , dict_value.len());
//                 if correct_count == 100 {
//                     // info!("Send dict_value Chunk to correct chanel: {:?}",dict_value);
//                     debug!(
//                         "Send dict_value Chunk reads {:?} to correct chanel: {:?}",
//                         correct_count, dict_key
//                     );
//                     dict_value.clear();
//                 }
//             }
//         }
//         // info!("readparts_idx: {:?}",readparts_idx);
//         // info!("readparts_dict_len: {:?}",readparts_dict.len());
//         // dbg!(record);
//     }
//     // 将所有readparts_dict内残留的未满的dict_value发送到correct channel
//     for (key, value) in readparts_dict.iter() {
//         if value.len() > 0 {
//             // info!("Send dict_value Chunk to correct chanel: {:?}",value);
//             debug!(
//                 "Last Send dict_value Chunk reads {:?} to correct chanel: {:?}",
//                 value.len(),
//                 key
//             );
//         }
//     }
// }
