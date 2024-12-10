use crate::bam::ReadPartChunk;
use crate::correcter;
use crate::correcter::Correcter;
use crate::poa::Aligner;
use crate::poa::POAGraph;
use anyhow::{Context, Result};
use bio::alignment::pairwise::Scoring;
use flume::Receiver;
use log::{debug, error, info};
use petgraph::dot::Dot;
use petgraph::graph;
use petgraph::{Directed, Graph, Incoming};
use std::thread;
use std::time::Instant;

pub struct BamHandler {
    pub correcter: Correcter,
    pub read_chunk: ReadPartChunk,
    pub graph: POAGraph,
    pub correct_ratio: f64,
    pub correct_fix_ratio: f64,
}
impl BamHandler {
    pub fn new(
        read_chunk: ReadPartChunk,
        correct_ratio: f64,
        correct_fix_ratio: f64,
    ) -> BamHandler {
        let mut read_chunk = read_chunk;
        let graph = Self::build_graph(&mut read_chunk);
        let correcter = Correcter::new("test".to_string(), graph.clone());
        BamHandler {
            correcter: correcter,
            read_chunk: read_chunk,
            graph: graph,
            correct_ratio: correct_ratio,
            correct_fix_ratio: correct_fix_ratio,
        }
    }

    pub fn build_graph(read_chunk: &mut ReadPartChunk) -> POAGraph {
        let first_read = &mut read_chunk[0];
        let scoring = Scoring::new(-1, 0, |a: u8, b: u8| if a == b { 3 } else { -4 });
        let mut aligner = Aligner::new(scoring, &first_read.read_seq);
        first_read.seq_node_index = aligner.first_seq_index.clone();
        // Process subsequent ReadInfo
        for read_info in &mut read_chunk[1..] {
            aligner
			.global(&read_info.read_seq);
        // if read_info.read_id == "TB2000CC3C-202406061734310_20240607005606_00517_000872770_11.62%4.2-F_3.7-R%00175.split.fastq.gz%-" {
            // 	
            // }
			read_info.seq_node_index = aligner.add_to_graph_and_correct();
            // if read_info.read_id == "TB2000CC3C-202406061734310_20240607020019_00729_003441578_12.17" {
                // if read_info.read_id == "TB2000CC3C-202406061734310_20240606183107_00177_000021521_11.46%4.2-F_3.7-R%00001.split.fastq.gz%-" {
                //     debug!("readid op: {:?}", aligner.alignment());
                //     debug!("readinfo: {:?}", read_info);
                //     debug!("readid seq: {:?}", std::str::from_utf8(&read_info.read_seq));
                // }
            //     info!("read_id:{} {:?}",read_info.read_id,aligner.alignment());
            // }
        }
        // info!("{:?}",read_chunk.ids());
        aligner.graph().clone()
    }
    pub fn correct_reads(&mut self) {
        for read_info in self.read_chunk.iter_mut() {
            // let (correct_indices, correct_qual) = self.correcter.correct_indices(
            //     &read_info.seq_node_index,
            //     &read_info.record.qual().to_owned(),
            // );
			let (correct_indices, correct_qual) = match self.correcter.correct_indices(
				&read_info.seq_node_index,
				&read_info.read_qual,
				self.correct_ratio,
				self.correct_fix_ratio,
			) {
				Ok((indices, qual)) => (indices, qual),
				Err(e) => {
					let g: Graph<char, i32, Directed, usize> = self.graph.map(|_, n| (*n) as char, |_, e| *e);
					debug!("Dot graph:\n{:?}", Dot::new(&g));
					// debug!(
					// 	"Graph: {:?}",
					// 	self.graph.map(|_, n| (*n) as char, |_, e| *e)
					// );
					// Optionally, log a portion of the data to avoid stack overflow
					// debug!(
					// 	"Error correcting indices for read_info with ID:{:?}\n{:?}\nNode Index: {:?} Seq: {:?}",
					// 	read_info.read_id,e,&read_info.seq_node_index, std::str::from_utf8(&read_info.read_seq)
					// );
					(
						read_info.seq_node_index.clone(),
						read_info.read_qual.clone(),
					)
				}
			};
            // let (correct_indices, correct_qual) = self.correcter.correct_indices(
            // 	&read_info.seq_node_index,
            // 	&read_info.read_qual.to_owned(),
            // 	self.correct_ratio,
            // 	self.correct_fix_ratio,
            // ).unwrap_or_else(|e| {
            // 	// Print debugging information when an error occurs
            // 	error!("Error correcting indices for read_info with ID: {}", read_info.read_id);
            // 	error!("Debugging Information:");
            // 	error!("Sequence Node Index: {:?}", read_info.seq_node_index);
            // 	error!("Quality Scores: {:?}", read_info.read_qual);
            // 	error!("Error: {:?}", e);
            // 	let g: Graph<char, i32, Directed, usize> = self.graph.map(|_, n| (*n) as char, |_, e| *e);
            // 	error!("Dot graph:\n{:?}", Dot::new(&g));
            // 	// Optionally, panic or handle the error in another way
            // 	error!("Failed to correct indices: {:?}", e);
            // 	(read_info.seq_node_index.to_owned(),read_info.read_qual.to_owned())
            // });

            // .expect(&format!("Error correcting indices for read_info with ID: {}", read_info.read_id));

            read_info.correct_seq = self.correcter.replace_with_node_weights(&correct_indices);
            read_info.correct_qual = correct_qual;
        }
    }
    pub fn info(&self) {
        let g: Graph<char, i32, Directed, usize> = self.graph.map(|_, n| (*n) as char, |_, e| *e);
        // debug!("Dot graph:\n{:?}", Dot::new(&g));
    }
}
pub fn correcter_receiver(
    rrx: Receiver<ReadPartChunk>,
    threads: usize,
    correct_ratio: f64,
    correct_fix_ratio: f64,
) -> Receiver<ReadPartChunk> {
    let (stx, srx) = flume::unbounded();
    for t in 0..threads {
        let start_time = Instant::now();
        let rrx = rrx.clone();
        let stx = stx.clone();
        thread::spawn(move || {
            let mut chunk_count = 0;
            let mut chunk_reads_num = 0;
            for readchunk in rrx.iter() {
                chunk_count += 1;
                // 过滤出需要校正的 reads
                let (correct_chunk, uncorrected_chunk): (Vec<_>, Vec<_>) =
                    readchunk.into_iter().partition(|rp| rp.need_correct());
                // 统计需要校正的 reads 数量
                // info!("correct_chunk.len() is {}; uncorrect_chunk.len() is {}", correct_chunk.len(), uncorrected_chunk.len());
                if !correct_chunk.is_empty() {
                    let correct_count = correct_chunk.len();
                    chunk_reads_num += correct_count;
                    let mut ch = BamHandler::new(correct_chunk, correct_ratio, correct_fix_ratio);
                    ch.correct_reads();
                    ch.info();
                    // info!("read1: {}", matched_reads.to_tsv());
                    stx.send(ch.read_chunk).expect("splitter send error");
                }
                if !uncorrected_chunk.is_empty() {
                    chunk_reads_num += uncorrected_chunk.len();
                    stx.send(uncorrected_chunk).expect("splitter send error");
                }
            }
            let elapsed_time = start_time.elapsed();
            info!(
                "threads {} process {} chunks {} reads. Time elapsed: {:.4?}",
                t, chunk_count, chunk_reads_num, elapsed_time
            )

        });
    }
    srx
}

#[test]
pub fn test_bam_handler() {}
