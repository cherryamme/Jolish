use crate::correcter;
use crate::correcter::Correcter;
use crate::fastq::ReadChunk;
use crate::poa::Aligner;
use crate::poa::POAGraph;
use bio::alignment::pairwise::Scoring;
use flume::Receiver;
use log::{debug, info};
use petgraph::dot::Dot;
use petgraph::graph;
use petgraph::{Directed, Graph, Incoming};
use std::thread;
use std::time::Instant;
use anyhow::{Context, Result};

pub struct CorrectHandler {
    pub correcter: Correcter,
    pub read_chunk: ReadChunk,
    pub graph: POAGraph,
    pub correct_ratio: f64,
    pub correct_fix_ratio: f64,
}
impl CorrectHandler {
    pub fn new(read_chunk: ReadChunk, correct_ratio: f64, correct_fix_ratio: f64) -> CorrectHandler {
        let mut read_chunk = read_chunk;
        let graph = Self::build_graph(&mut read_chunk);
        let correcter = Correcter::new("test".to_string(), graph.clone());
        CorrectHandler {
            correcter: correcter,
            read_chunk: read_chunk,
            graph: graph,
            correct_ratio: correct_ratio,
            correct_fix_ratio: correct_fix_ratio,
        }
    }

    pub fn build_graph(read_chunk: &mut ReadChunk) -> POAGraph {
        let first_read = &mut read_chunk.reads[0];
        let scoring = Scoring::new(-1, 0, |a: u8, b: u8| if a == b { 3 } else { -4 });
        let mut aligner = Aligner::new(scoring, first_read.record.seq());
        first_read.seq_node_index = aligner.first_seq_index.clone();
        // Process subsequent ReadInfo
        for read_info in &mut read_chunk.reads[1..] {
            read_info.seq_node_index = aligner
                .global(read_info.record.seq())
                .add_to_graph_and_correct();
            // if read_info.read_id == "TB2000CC3C-202406061734310_20240607020019_00729_003441578_12.17" {
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
            let (correct_indices, correct_qual) = self.correcter.correct_indices(
                &read_info.seq_node_index,
                &read_info.record.qual().to_owned(),
                self.correct_ratio,
                self.correct_fix_ratio,
            ).unwrap_or_else(|e| {
                // Print debugging information when an error occurs
                println!("Error correcting indices for read_info with ID: {}", read_info.read_id);
                println!("Debugging Information:");
                println!("Sequence Node Index: {:?}", read_info.seq_node_index);
                println!("Quality Scores: {:?}", read_info.record.qual());
                println!("Error: {:?}", e);
                let g: Graph<char, i32, Directed, usize> = self.graph.map(|_, n| (*n) as char, |_, e| *e);
                info!("Dot graph:\n{:?}", Dot::new(&g));
                // Optionally, panic or handle the error in another way
                panic!("Failed to correct indices: {:?}", e);
            });
            
            // .expect(&format!("Error correcting indices for read_info with ID: {}", read_info.read_id));

            read_info.correct_seq = self.correcter.replace_with_node_weights(&correct_indices);
            read_info.correct_qual = correct_qual;
        }
    }
    pub fn info(&self) {
        let g: Graph<char, i32, Directed, usize> = self.graph.map(|_, n| (*n) as char, |_, e| *e);
        debug!("Dot graph:\n{:?}", Dot::new(&g));
    }
}

pub fn correcter_receiver(
    rrx: Receiver<ReadChunk>,
    threads: usize,
    chunksize: usize,
    correct_ratio: f64,
    correct_fix_ratio: f64,
) -> Receiver<ReadChunk> {
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
                chunk_reads_num += readchunk.reads_num();
                let mut ch = CorrectHandler::new(readchunk, correct_ratio, correct_fix_ratio);
                ch.correct_reads();
                ch.info();
                // info!("read1: {}", matched_reads.to_tsv());
                stx.send(ch.read_chunk).expect("splitter send error");
            }
            let elapsed_time = start_time.elapsed();
            info!(
                "threads {} process {} chunks {} reads. Time elapsed: {:.4?}",
                t,
                chunk_count,
                chunk_reads_num,
                elapsed_time
            )
        });
    }
    srx
}
