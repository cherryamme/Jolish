use bio::io::fastq::{Reader, Record};
use flate2::read::MultiGzDecoder;
use flume::{unbounded, Sender, Receiver};
use log::{info,debug,error};
use petgraph::graph::NodeIndex;
use std::clone;
use std::ffi::OsStr;
use std::{
    fs::File,
    io::{BufReader, Read},
    path::PathBuf,
};
use std::time::Instant;
use std::collections::HashMap;

pub struct ReadInfo {
    pub record: Record,
    pub strand_orient: String,
    pub read_id: String,
    pub read_type: String,
    pub write_to_fq: bool,
    pub seq_node_index: Vec<NodeIndex<usize>>,
    pub correct_seq: Vec<u8>,
    pub correct_qual: Vec<u8>,
    pub read_len: usize,
}
impl ReadInfo {
    pub fn new(record: Record, split_param: char, type_index: usize, orient_index: usize ) -> ReadInfo {
        let mut readinfo = ReadInfo {
            record: record.clone(),
            read_id: String::from("unknown"),
            strand_orient: String::from("unknown"),
            read_type: String::from("default"),
            write_to_fq: false,
            seq_node_index: Vec::new(),
            correct_seq: Vec::new(),
            correct_qual: Vec::new(),
            read_len: record.seq().len(),
        };
        readinfo.parser_id(split_param, type_index, orient_index);
        readinfo
    }
    fn parser_id(&mut self, split_param: char, type_index: usize, orient_index: usize) {
        let id_parts: Vec<&str> = self.record.id().split(split_param).collect();
        self.read_id = id_parts[0].to_string();
        let type_part = id_parts.get(type_index).copied();
        let orient_part = id_parts.get(orient_index).copied();
        if let Some(type_part) = type_part {
            self.read_type = type_part.to_string();
        }
        if let Some(orient_part) = orient_part {
            self.strand_orient = orient_part.to_string();
        }
    }
    pub fn info(&self) -> String {
        self.record.id().to_string()
    }
}

pub struct ReadChunk {
    pub reads: Vec<ReadInfo>,
    key: (String, String),
}

impl ReadChunk {
    fn new(capacity: usize, key: (String, String)) -> Self {
        Self {
            reads: Vec::with_capacity(capacity),
            key,
        }
    }
    // 添加readinfo，并返回是否已满
    fn add(&mut self, read_info: ReadInfo) -> bool {
        self.reads.push(read_info);
        self.reads.len() == self.reads.capacity()
    }
    pub fn info(&self) -> (String, String) {
        self.key.clone()
    }
    // 返回可变引用的迭代器
    pub fn iter_mut(&mut self) -> std::slice::IterMut<ReadInfo> {
        self.reads.iter_mut()
    }
    // 返回所有reads的id名称
    pub fn ids(&self) -> Vec<String> {
        self.reads.iter().map(|read| read.read_id.clone()).collect()
    }
    pub fn reads_num(&self) -> usize {
        self.reads.len()
    }
}



const BUFSIZE: usize = 10 * 1024 * 1024;

fn is_gz(path: &PathBuf) -> bool {
    match path.extension().and_then(OsStr::to_str) {
        Some(ext) => ext == "gz",
        None => false,
    }
}

pub fn spawn_reader(files: Vec<String>, chunksize: usize, split_param: char, type_index: usize, orient_index: usize, max_reads: Option<usize>) -> Receiver<ReadChunk> {
    let (rtx, rrx) = unbounded();
    std::thread::spawn(move || {
        let start_time = Instant::now();
        if files.is_empty() {
            info!("no input file, loading from stdin...");
            let stdin_handle = std::io::stdin();
            process_file(stdin_handle, &rtx, None, chunksize, split_param, type_index, orient_index, max_reads);
        } else {
            for file in files {
                let path = PathBuf::from(&file);
                if path.exists() {
                    let raw_handle = File::open(&path)
                        .expect(format!("Error opening input: {}", path.display()).as_str());
                    process_file(raw_handle, &rtx, Some(path), chunksize, split_param, type_index, orient_index, max_reads);
                } else {
                    panic!("File {} does not exist", path.display());
                }
            }
        }


        let elapsed_time = start_time.elapsed();
        info!("Loading Reads data done! Time elapsed: {:.4?}", elapsed_time)
    });
    rrx
}


fn process_file<R: Read + 'static>(handle: R, rtx: &Sender<ReadChunk>, path: Option<PathBuf>, chunksize: usize, split_param: char, type_index: usize, orient_index: usize, max_reads: Option<usize>) {
    // Process a file handle, either from a file or stdin
    // If the file is a gzip file, use MultiGzDecoder to decode it
    // If the file is a fastq file, use BufReader to read it
    // If the file is stdin, use BufReader to read it
    // Send the read info to the channel
    let buf_handle = BufReader::with_capacity(BUFSIZE, handle);
    let maybe_decoder_handle = {
        if let Some(path) = path {
            if is_gz(&path) {
                info!("loading gzip file:{:?}", path);
                Box::new(MultiGzDecoder::new(buf_handle)) as Box<dyn Read>
            } else {
                info!("loading fastq file:{:?}", path);
                Box::new(buf_handle) as Box<dyn Read>
            }
        } else {
            Box::new(buf_handle) as Box<dyn Read>
        }
    };
    let fastq_reader = Reader::new(maybe_decoder_handle);
    let mut chunks: HashMap<(String, String), ReadChunk> = HashMap::new();
    // 将发送单条reads修改为发送ReadChunk
    // for record in fastq_reader.records() {
    //     let readinfo = ReadInfo::new(record.unwrap());
    //     rtx.send(readinfo).expect("Error sending");
    // }
    let max_reads_limit = max_reads.unwrap_or(usize::MAX);
    let mut process_num = 0;
    for record in fastq_reader.records() {
        if process_num >= max_reads_limit {break};
        process_num += 1;
        let record = match record {
            Ok(record) => record,
            Err(e) => {
                error!("Error: {}", e);
                continue;
            }
        };    
        // debug!("record id is {}",record.clone().id());
        let read_info = ReadInfo::new(record, split_param, type_index, orient_index);
        let key = (read_info.read_type.clone(), read_info.strand_orient.clone());
        let chunk = chunks.entry(key.clone()).or_insert_with(|| ReadChunk::new(chunksize, key.clone()));
        if chunk.add(read_info) {
            let key = key.clone();
            let chunk_to_send = chunks.remove(&key).unwrap();
            rtx.send(chunk_to_send).expect("Error sending");
        }

    }
    // Send any remaining reads in the last chunk for each group
    for (_, read_chunk) in chunks {
        if !read_chunk.reads.is_empty() {
            rtx.send(read_chunk).expect("Error sending");
        }
    }
}



#[test]
fn test_fastq() {
    pretty_env_logger::init();
    let mut files = Vec::new();
    files.push("test.fq".to_string());
    let rrx: Receiver<ReadChunk> = spawn_reader(files, 10, '%', 1, 3, Some(200));
    let mut count = 0;
    for mut chunk in rrx.iter() {
        count += 1;
        debug!("No {} chunk length is {}",count,chunk.reads.len());
        debug!("key is {:?}",chunk.info());
        for read in chunk.iter_mut() {
            debug!("read id is {}",read.info());
        }
    }
    info!("chunk length is {}",count);
}
