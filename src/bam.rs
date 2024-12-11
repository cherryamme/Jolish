use flate2::read;
use flume::{unbounded, Receiver, Sender};
use log::{debug, error, info};
use petgraph::graph::NodeIndex;
use rust_htslib::bam::record;
use rust_htslib::bam::{self, ext::BamRecordExtensions, Read, Record};
use rust_htslib::bcf::record::Info;
use serde::de::value;
use std::collections::HashMap;
use std::sync::{Arc, Mutex};
use std::time::Instant;
use std::{clone, io};
use std::{fmt, str::FromStr};
use std::{fs::File, io::BufReader, path::PathBuf};

pub type ReadPartChunk = Vec<ReadPart>;

#[derive(Debug, Clone)]
struct Region {
    chrom: String,
    start: i32,
    end: i32,
}
impl Region {
    pub fn new(chrom: String, start: i32, end: i32) -> Self {
        Self { chrom, start, end }
    }
}
impl std::str::FromStr for Region {
    type Err = std::num::ParseIntError;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let parts: Vec<&str> = s.split(':').collect();
        let chrom = parts[0].to_owned();
        let parts: Vec<&str> = parts[1].split('-').collect();
        let start = parts[0].parse()?;
        let end = parts[1].parse()?;
        Ok(Region { chrom, start, end })
    }
}
#[derive(Debug)]
pub struct ReadpartsNumDictValue {
    pub orient: String,
    pub record_info: Vec<RecordInfo>,
    pub record_num: usize,
    pub record_count: usize,
    pub chunk_num: usize,
}

#[derive(Debug, Clone)]
pub struct RecordInfo {
    pub record_id: String,
    pub record_order: usize,
    pub orient: String,
}

#[derive(Debug)]
pub struct ReadParts_idx {
    pub chrom: usize,
    pub windows_num: i64,
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
        read_qual: &[u8],
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

pub fn get_readparts_idx(record: &Record, chunk_window: i64) -> Vec<ReadParts_idx> {
    let infos = record.aligned_pairs();
    let softclip_len = record.cigar().leading_softclips();
    let mut chunk_start = softclip_len;
    let mut windows_num = record.pos() / chunk_window;
    let mut reader_order = (softclip_len / chunk_window) as usize;
    let mut readparts_idxs = Vec::new();
    // debug!("leading_softclips:{:?}  record.pos: {:?} reference_start:{:?}",record.cigar().leading_softclips(),record.pos(),record.reference_start());
    for info in infos {
        
        let read_pos = info[0];
        let ref_pos = info[1];
        if ref_pos / chunk_window > windows_num {
            let readparts_idx = ReadParts_idx {
                chrom: record.tid() as usize,
                windows_num,
                read_order: reader_order,
                start: chunk_start as usize,
                end: read_pos as usize,
                ref_pos: ref_pos as usize,
            };
            readparts_idxs.push(readparts_idx);
            windows_num = ref_pos / chunk_window;
            reader_order += 1;
            chunk_start = read_pos;
        }
    }
    // debug!("readparts_idxs: {:?}",readparts_idxs);
    readparts_idxs
}





fn parse_satag_to_count(record: &Record) -> usize {
    match record.aux(b"SA"){
        Ok(value) => {
            let content = format!("{:?}", value);
            content.split(";").count()
        }
        Err(e) => {
            1
        }
    }
}













pub fn spawn_bam_reader(
    files: Vec<String>,
    chunk_size: usize,
    chunk_window: usize,
    max_reads: Option<usize>,
    region: Option<String>,
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
            // log日志，打印bam文件的header信息
            if let Some(region) = region.clone() {
                // 解析需要提取的region
                let mut bam = bam::IndexedReader::from_path(&file).unwrap();
                let header = bam::Header::from_template(bam.header());
                info!("Process bam {:?} in region: {:?}", file, region);
                let region = Region::from_str(&region).unwrap();
                // 获取bam文件的染色体信息; 根据需要进行片段截取
                bam.fetch((&region.chrom, region.start, region.end))
                    .unwrap();
            }
            let mut bam = bam::Reader::from_path(&file).unwrap();
            let header = bam::Header::from_template(bam.header());
            // 循环处理bam文件中的每一行
            info!("Process bam: {:?}", file);
            for result in bam.records() {
                if let Some(max) = max_reads {
                    if record_count >= max {
                        break;
                    }
                }
                // 处理的record数量加1
                let record = result.unwrap();
                // 如果record是unmapped的，跳过
                if record.is_unmapped() {
                    continue;
                }
                record_count += 1;
                // 读取record是正向还是反向
                let orient = if record.is_reverse() {
                    "reverse"
                } else {
                    "forward"
                };
                // info!("record_SA: {:?}", parse_satag_to_count(&record));
                let sa_tag_count = parse_satag_to_count(&record);
                // TODO sofclip不支持校正，需要修改
                // 如果比对上多个地方，直接跳过，进行下一个record的校正
                // if sa_tag_count > 1 {
                //     continue;
                // }
                // 计算record的分窗信息
                let readparts_idxs = get_readparts_idx(&record,100);
                // 解析record的信息
                let read_id: String = std::str::from_utf8(record.qname()).unwrap().to_owned();
                let read_seq = &record.seq().as_bytes();
                let read_qual: &[u8] = record.qual();
                // 生成record_info结构体，内部包含record_id, record_order, orient
                let record_info = RecordInfo {
                    record_id: read_id.clone(),
                    // TODO record_order用于标记record的顺序，用于处理多个比对的read
                    record_order: 0,
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
                        value.chunk_num += readparts_idxs.len(); // 假设每次调用增加一个块
                    })
                    .or_insert(ReadpartsNumDictValue {
                        record_info: vec![record_info],
                        orient: orient.to_string(),
                        record_num: sa_tag_count,
                        record_count: 1,                // 初始化为1
                        chunk_num: readparts_idxs.len(), // 初始化为1
                    });
                    // 根据record的分窗信息，生成readpart，并将readpart发送到rtx
                    for readpart_idx in readparts_idxs.iter() {
                        // debug!("readpart_idx: {:?}",readpart_idx);
                        // 分窗内创建readpart结构
                        let readpart =
                            ReadPart::new(read_id.clone(),read_id.clone(), read_seq, read_qual, readpart_idx);

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

// #[test]
// pub fn test_spawn_bam_reader() {
//     pretty_env_logger::init();
//     let files = vec!["/home/jiangchen/project/Jolish/example/bam_example.fq.bam".to_string()];
//     let rrx = spawn_bam_reader(files, 100, 100, Some(200));
//     // for readparts in rrx.iter() {
//     //     info!("Readparts: {:?}",readparts);
//     // }
// }

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
