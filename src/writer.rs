// use std::collections::HashMap;
use std::fs::File;
use std::io::Write;
use flate2::write::GzEncoder;
use flate2::Compression;
use log::info;
// use std::io::Result;
use std::path::Path;
// use bio::io::fastq::Writer;
use std::fs::create_dir_all;
use crate::fastq::ReadInfo;
use crate::fastq::ReadChunk;
use std::io::BufWriter;
// use std::thread;
use flume::{Receiver, Sender, unbounded};
use std::time::Instant;
use std::thread;


// pub struct WriterManager {
//     sender: Option<Sender<ReadInfo>>,
//     handle: Option<thread::JoinHandle<()>>,
// }

// impl WriterManager {
//     pub fn new(outfile: String) -> WriterManager {
//         info!("Creating writer manager, start writing...");
//         let (tx, rx) = unbounded();
//         let filepath_string = format!("{}.fq.gz", outfile);
//         let filepath = Path::new(&filepath_string);
//         let filedir = &filepath.parent().unwrap();
//         create_dir_all(&filedir).expect("fail to create output directory");
//         let file = File::create(&filepath).expect("fail to create output fq.gz");
//         let encoder = GzEncoder::new(file, Compression::default());
//         let writer = BufWriter::with_capacity(1_000_000, encoder);
//         let handle = Self::start_writing_thread(writer, rx);

//         WriterManager {
//             sender: Some(tx),
//             handle: Some(handle),
//         }
//     }

//     pub fn write(&mut self, readinfo: ReadInfo) -> Result<()> {
//         if !readinfo.write_to_fq {
//             return Ok(());
//         }
//         if let Some(sender) = &self.sender {
//             sender.send(readinfo).expect("readinfo to writer send fail");
//         }
//         Ok(())
//     }

//     fn start_writing_thread(mut writer: BufWriter<GzEncoder<File>>, rx: Receiver<ReadInfo>) -> thread::JoinHandle<()> {
//         thread::spawn(move || {
//             for readinfo in rx.iter() {
//                 let id = readinfo.out_record.id();
//                 let seq = std::str::from_utf8(readinfo.out_record.seq()).expect("Not a valid UTF-8 sequence");
//                 let qual = std::str::from_utf8(readinfo.out_record.qual()).expect("Not a valid UTF-8 sequence");
//                 let record_str = format!("@{}\n{}\n+\n{}\n", id, seq, qual);
//                 write!(writer, "{}", record_str).unwrap();
//             }
//         })
//     }

//     pub fn drop(&mut self) {
//         // When the `Sender`s are dropped, the corresponding writing threads will receive a `Disconnected` error and exit.
//         info!("Writing fastq.gz. May cost some time..");
//         self.writers.clear();
//         // Wait for all writing threads to finish.
//         for handle in self.handles.drain(..) {
//             handle.join().expect("Writing thread panicked");
//         }
//     }
// }



// pub fn write_to_fq(filepath: &str, readinfo: ReadInfo){
//     let (stx, srx) = flume::unbounded();
//     let filepath = Path::new(filepath);
//     let filedir = &filepath.parent().unwrap();
//     create_dir_all(&filedir).expect("fail to create output directory");
//     let file = File::create(&filepath).expect("fail to create output fq.gz");
//     let encoder = GzEncoder::new(file, Compression::default());
//     let mut writer = BufWriter::with_capacity(1_000_000, encoder);
//     let id =  readinfo.record.id();
//     let seq =  std::str::from_utf8(&readinfo.correct_seq).expect("Not a valid UTF-8 sequence");
//     let qual = std::str::from_utf8(&readinfo.correct_qual).expect("Not a valid UTF-8 sequence");
//     let record_str = format!("@{}\n{}\n+\n{}\n", id, seq, qual);
//     write!(writer, "{}", record_str).unwrap();
// }


pub fn writer_receiver(rrx: Receiver<ReadChunk>, outfile: &str){
    let start_time = Instant::now();
    let filepath = Path::new(outfile);
    let filedir = &filepath.parent().unwrap();
    create_dir_all(&filedir).expect("fail to create output directory");
    let file = File::create(&filepath).expect("fail to create output fq.gz");
    let encoder = GzEncoder::new(file, Compression::default());
    let mut writer = BufWriter::with_capacity(1_000_000, encoder);
    for readchunk in rrx{
            for readinfo in readchunk.reads{
                let id =  readinfo.record.id();
                let seq =  std::str::from_utf8(&readinfo.correct_seq).expect("Not a valid UTF-8 sequence");
                let qual = std::str::from_utf8(&readinfo.correct_qual).expect("Not a valid UTF-8 sequence");
                let record_str = format!("@{}\n{}\n+\n{}\n", id, seq, qual);
                write!(writer, "{}", record_str).unwrap();
            }
        }
    let elapsed_time = start_time.elapsed();
    info!("writer to fastq {}.\n Time elapsed: {:.4?}",outfile,elapsed_time)
}
