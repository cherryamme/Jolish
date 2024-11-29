
use bio::alignment::pairwise::{MatchFunc, Scoring};
use bio::bio_types::sequence;
use bio::utils::TextSlice;
use petgraph::graph::{self, NodeIndex};
use petgraph::visit::Topo;
use petgraph::dot::Dot;
use petgraph::{Directed, Graph, Incoming};
use std::cmp::{max, Ordering};
use serde::{Serialize, Deserialize};
use std::str::from_utf8;
pub type POAGraph = Graph<u8, i32, Directed, usize>;
pub const MIN_SCORE: i32 = -858_993_459; // negative infinity; see alignment/pairwise/mod.rs

#[derive(Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Serialize, Deserialize)]
pub enum AlignmentOperation {
    Match(Option<(usize, usize, char)>), // added char for base type
    Del(Option<(usize, usize, char)>),   // added char for base type
    Ins(Option<(usize, char)>),          // added char for base type
    Xclip(usize),
    Yclip(usize, usize), // to, from
}


#[derive(Copy, Clone, Debug, Serialize, Deserialize)]
pub struct TracebackCell {
    score: i32,
    op: AlignmentOperation,
}
impl Ord for TracebackCell {
    fn cmp(&self, other: &TracebackCell) -> Ordering {
        self.score.cmp(&other.score)
    }
}

impl PartialOrd for TracebackCell {
    fn partial_cmp(&self, other: &TracebackCell) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl PartialEq for TracebackCell {
    fn eq(&self, other: &TracebackCell) -> bool {
        self.score == other.score
    }
}

//impl Default for TracebackCell { }

impl Eq for TracebackCell {}

#[derive(Default, Clone, Eq, PartialEq, Ord, PartialOrd, Debug)]
pub struct Traceback {
    rows: usize,
    cols: usize,

    // store the last visited node in topological order so that
    // we can index into the end of the alignment when we backtrack
    last: NodeIndex<usize>,
    matrix: Vec<(Vec<TracebackCell>, usize, usize)>,
}

impl Traceback {
	fn new() -> Self {
        Traceback {
            rows: 0,
            cols: 0,
            last: NodeIndex::new(0),
            matrix: Vec::new(),
        }
    }
	fn with_capacity(m: usize, n: usize) -> Self {
        // each row of matrix contain start end position and vec of traceback cells
        let matrix: Vec<(Vec<TracebackCell>, usize, usize)> = vec![(vec![], 0, n + 1); m + 1];
        Traceback {
            rows: m,
            cols: n,
            last: NodeIndex::new(0),
            matrix,
        }
    }
	fn initialize_scores(&mut self, gap_open: i32, yclip: i32) {
        for j in 0..=self.cols {
            self.matrix[0].0.push(max(
                TracebackCell {
                    score: (j as i32) * gap_open,
                    op: AlignmentOperation::Ins(None),
                },
                TracebackCell {
                    score: yclip,
                    op: AlignmentOperation::Yclip(0, j),
                },
            ));
        }
        self.matrix[0].0[0] = TracebackCell {
            score: 0,
            op: AlignmentOperation::Match(None),
        };
    }

	fn set(&mut self, i: usize, j: usize, cell: TracebackCell) {
        // set the matrix cell if in band range
        if !(self.matrix[i].1 > j || self.matrix[i].2 < j) {
            let real_position = j - self.matrix[i].1;
            self.matrix[i].0[real_position] = cell;
        }
    }
    fn get(&self, i: usize, j: usize) -> &TracebackCell {
        // get the matrix cell if in band range else return the appropriate values
        if !(self.matrix[i].1 > j || self.matrix[i].2 <= j || self.matrix[i].0.is_empty()) {
            let real_position = j - self.matrix[i].1;
            &self.matrix[i].0[real_position]
        }
        // behind the band, met the edge
        else if j == 0 {
            &TracebackCell {
                score: MIN_SCORE,
                op: AlignmentOperation::Del(None),
            }
        }
        // infront of the band
        else if j >= self.matrix[i].2 {
            &TracebackCell {
                score: MIN_SCORE,
                op: AlignmentOperation::Ins(None),
            }
        }
        // behind the band
        else {
            &TracebackCell {
                score: MIN_SCORE,
                op: AlignmentOperation::Match(None),
            }
        }
    }
	fn new_row(
        &mut self,
        row: usize,
        size: usize,
        gap_open: i32,
        xclip: i32,
        start: usize,
        end: usize,
    ) {
        self.matrix[row].1 = start;
        self.matrix[row].2 = end;
        // when the row starts from the edge
        if start == 0 {
            self.matrix[row].0.push(max(
                TracebackCell {
                    score: (row as i32) * gap_open,
                    op: AlignmentOperation::Del(None),
                },
                TracebackCell {
                    score: xclip,
                    op: AlignmentOperation::Xclip(0),
                },
            ));
        } else {
            self.matrix[row].0.push(TracebackCell {
                score: MIN_SCORE,
                op: AlignmentOperation::Match(None),
            });
        }
        for _ in 1..=size {
            self.matrix[row].0.push(TracebackCell {
                score: MIN_SCORE,
                op: AlignmentOperation::Match(None),
            });
        }
    }

	pub fn alignment(&self) -> Alignment {
        // optimal AlignmentOperation path
        let mut ops: Vec<AlignmentOperation> = vec![];

        // Now backtrack through the matrix to construct an optimal path
        let mut i = self.last.index() + 1;
        let mut j = self.cols;

        while i > 0 || j > 0 {
            // push operation and edge corresponding to (one of the) optimal
            // routes
            ops.push(self.get(i, j).op);
            match self.get(i, j).op {
                AlignmentOperation::Match(Some((p, _, _))) => {
                    i = p + 1;
                    j -= 1;
                }
                AlignmentOperation::Del(Some((p, _, _))) => {
                    i = p + 1;
                }
                AlignmentOperation::Ins(Some((p, _))) => {
                    i = p + 1;
                    j -= 1;
                }
                AlignmentOperation::Match(None) => {
                    i = 0;
                    j -= 1;
                }
                AlignmentOperation::Del(None) => {
                    i -= 1;
                }
                AlignmentOperation::Ins(None) => {
                    j -= 1;
                }
                AlignmentOperation::Xclip(r) => {
                    i = r;
                }
                AlignmentOperation::Yclip(r, _) => {
                    j = r;
                }
            }
        }

        ops.reverse();

        Alignment {
            score: self.get(self.last.index() + 1, self.cols).score,
            operations: ops,
        }
    }
}

#[derive(Default, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Serialize, Deserialize)]
pub struct Alignment {
    pub score: i32,
    //    xstart: Edge,
    operations: Vec<AlignmentOperation>,
}

#[derive(Default, Clone, Debug)]
pub struct Poa<F: MatchFunc> {
    scoring: Scoring<F>,
    pub graph: POAGraph,
}
impl<F: MatchFunc> Poa<F> {
	pub fn from_string(scoring: Scoring<F>, seq: TextSlice) -> Self {
        let mut graph: Graph<u8, i32, Directed, usize> =
		Graph::with_capacity(seq.len(), seq.len() - 1);
        let mut prev: NodeIndex<usize> = graph.add_node(seq[0]);
        let mut node: NodeIndex<usize>;
        for base in seq.iter().skip(1) {
            node = graph.add_node(*base);
            graph.add_edge(prev, node, 1);
            prev = node;
        }
        Poa { scoring, graph }
    }
	pub fn from_string_get_index(scoring: Scoring<F>, seq: TextSlice) -> (Self, Vec<NodeIndex<usize>>) {
        let mut graph: Graph<u8, i32, Directed, usize> =
            Graph::with_capacity(seq.len(), seq.len() - 1);
        let mut prev: NodeIndex<usize> = graph.add_node(seq[0]);
        let mut node: NodeIndex<usize>;
        let mut node_indices: Vec<NodeIndex<usize>> = vec![prev]; // 初始化第一个节点索引
        
        for base in seq.iter().skip(1) {
            node = graph.add_node(*base);
            graph.add_edge(prev, node, 1);
            node_indices.push(node); // 添加节点索引到向量
            prev = node;
        }
        
        (Poa { scoring, graph }, node_indices)
    }
	pub fn custom(&self, query: TextSlice) -> Traceback {
		assert!(self.graph.node_count() != 0);
		// dimensions of the traceback matrix
		let (m, n) = (self.graph.node_count(), query.len());
		// save score location of the max scoring node for the query for suffix clipping
		let mut max_in_column = vec![(0, 0); n + 1];
		let mut traceback = Traceback::with_capacity(m, n);
		traceback.initialize_scores(self.scoring.gap_open, self.scoring.yclip_prefix);
		// construct the score matrix (O(n^2) space)
		let mut topo = Topo::new(&self.graph);
		while let Some(node) = topo.next(&self.graph) {
			// reference base and index
			let r = self.graph.raw_nodes()[node.index()].weight; // reference base at previous index
			let r_char = char::from(r);
			let i = node.index() + 1; // 0 index is for initialization so we start at 1
			traceback.last = node;
			// iterate over the predecessors of this node
			let prevs: Vec<NodeIndex<usize>> =
				self.graph.neighbors_directed(node, Incoming).collect();
			traceback.new_row(
				i,
				n + 1,
				self.scoring.gap_open,
				self.scoring.xclip_prefix,
				0,
				n + 1,
			);
			// query base and its index in the DAG (traceback matrix rows)
			for (query_index, query_base) in query.iter().enumerate() {
				let j = query_index + 1; // 0 index is initialized so we start at 1
										 // match and deletion scores for the first reference base
				let max_cell = if prevs.is_empty() {
					TracebackCell {
						score: traceback.get(0, j - 1).score
							+ self.scoring.match_fn.score(r, *query_base),
						op: AlignmentOperation::Match(None),
					}
				} else {
					let mut max_cell = max(
						TracebackCell {
							score: MIN_SCORE,
							op: AlignmentOperation::Match(None),
						},
						TracebackCell {
							score: self.scoring.xclip_prefix,
							op: AlignmentOperation::Xclip(0),
						},
					);
					for prev_node in &prevs {
						let i_p: usize = prev_node.index() + 1; // index of previous node
						max_cell = max(
							max_cell,
							max(
								TracebackCell {
									score: traceback.get(i_p, j - 1).score
										+ self.scoring.match_fn.score(r, *query_base),
									op: AlignmentOperation::Match(Some((i_p - 1, i - 1, r_char))),
								},
								TracebackCell {
									score: traceback.get(i_p, j).score + self.scoring.gap_open,
									op: AlignmentOperation::Del(Some((i_p - 1, i, r_char))),
								},
							),
						);
					}
					max_cell
				};
				let score = max(
					max_cell,
					TracebackCell {
						score: traceback.get(i, j - 1).score + self.scoring.gap_open,
						op: AlignmentOperation::Ins(Some((i - 1, r_char))),
					},
				);
				traceback.set(i, j, score);
				if max_in_column[j].0 < score.score {
					max_in_column[j].0 = score.score;
					max_in_column[j].1 = i;
				}
			}
		}
		// X suffix clipping
		let mut max_in_row = (0, 0);
		for (col_index, &(score, col_max_row)) in max_in_column.iter().enumerate() {
			// avoid pointing to itself
			if col_max_row == traceback.last.index() + 1 {
				continue;
			}
			let maxcell = max(
				*traceback.get(traceback.last.index() + 1, col_index),
				TracebackCell {
					score: score + self.scoring.xclip_suffix,
					op: AlignmentOperation::Xclip(col_max_row),
				},
			);
			if max_in_row.0 < maxcell.score {
				max_in_row.0 = maxcell.score;
				max_in_row.1 = col_index;
			}
			traceback.set(traceback.last.index() + 1, col_index, maxcell);
		}
		// Y suffix clipping from the last node
		let maxcell = max(
			*traceback.get(traceback.last.index() + 1, n),
			TracebackCell {
				score: max_in_row.0 + self.scoring.yclip_suffix,
				op: AlignmentOperation::Yclip(max_in_row.1, n),
			},
		);
		if max_in_row.1 != n {
			traceback.set(traceback.last.index() + 1, n, maxcell);
		}
		traceback
	}
	pub fn add_alignment(&mut self, aln: &Alignment, seq: TextSlice) {
        let head = Topo::new(&self.graph).next(&self.graph).unwrap();
        let mut prev: NodeIndex<usize> = NodeIndex::new(head.index());
        let mut i: usize = 0;
        let mut edge_not_connected: bool = false;
        for op in aln.operations.iter() {
            match op {
                AlignmentOperation::Match(None) => {
                    let node: NodeIndex<usize> = NodeIndex::new(head.index());
                    if (seq[i] != self.graph.raw_nodes()[head.index()].weight) && (seq[i] != b'X') {
                        let node = self.graph.add_node(seq[i]);
                        if edge_not_connected {
                            self.graph.add_edge(prev, node, 1);
                        }
                        edge_not_connected = false;
                        prev = node;
                    }
                    if edge_not_connected {
                        self.graph.add_edge(prev, node, 1);
                        prev = node;
                        edge_not_connected = false;
                    }
                    i += 1;
                }
                AlignmentOperation::Match(Some((_, p, _))) => {
                    let node = NodeIndex::new(*p);
                    if (seq[i] != self.graph.raw_nodes()[*p].weight) && (seq[i] != b'X') {
                        let node = self.graph.add_node(seq[i]);
                        self.graph.add_edge(prev, node, 1);
                        prev = node;
                    } else {
                        // increment node weight
                        match self.graph.find_edge(prev, node) {
                            Some(edge) => {
                                *self.graph.edge_weight_mut(edge).unwrap() += 1;
                            }
                            None => {
                                if prev.index() != head.index() && prev.index() != node.index() {
                                    self.graph.add_edge(prev, node, 1);
                                }
                            }
                        }
                        prev = NodeIndex::new(*p);
                    }
                    i += 1;
                }
                AlignmentOperation::Ins(None) => {
                    let node = self.graph.add_node(seq[i]);
                    if edge_not_connected {
                        self.graph.add_edge(prev, node, 1);
                    }
                    prev = node;
                    edge_not_connected = true;
                    i += 1;
                }
                AlignmentOperation::Ins(Some(_)) => {
                    let node = self.graph.add_node(seq[i]);
                    self.graph.add_edge(prev, node, 1);
                    prev = node;
                    i += 1;
                }
                AlignmentOperation::Del(_) => {} // we should only have to skip over deleted nodes and xclip
                AlignmentOperation::Xclip(_) => {}
                AlignmentOperation::Yclip(_, r) => {
                    i = *r;
                }
            }
        }
    }
	pub fn add_alignment_get_index(&mut self, aln: &Alignment, seq: TextSlice) -> Vec<NodeIndex<usize>> {
		let head = Topo::new(&self.graph).next(&self.graph).unwrap();
		let mut prev: NodeIndex<usize> = NodeIndex::new(head.index());
		let mut i: usize = 0;
		let mut edge_not_connected: bool = false;
		let mut node_indices: Vec<NodeIndex<usize>> = Vec::new(); // 新增的向量，用于存储节点索引
	
		for op in aln.operations.iter() {
			match op {
				AlignmentOperation::Match(None) => {
					let node: NodeIndex<usize> = NodeIndex::new(head.index());
					if (seq[i] != self.graph.raw_nodes()[head.index()].weight) && (seq[i] != b'X') {
						let node = self.graph.add_node(seq[i]);
						if edge_not_connected {
							self.graph.add_edge(prev, node, 1);
						}
						edge_not_connected = false;
						prev = node;
					}
					if edge_not_connected {
						self.graph.add_edge(prev, node, 1);
						prev = node;
						edge_not_connected = false;
					}
					node_indices.push(node); // 添加节点索引到向量
					i += 1;
				}
				AlignmentOperation::Match(Some((_, p, _))) => {
					let node = NodeIndex::new(*p);
					if (seq[i] != self.graph.raw_nodes()[*p].weight) && (seq[i] != b'X') {
						let node = self.graph.add_node(seq[i]);
						self.graph.add_edge(prev, node, 1);
						prev = node;
					} else {
						// increment node weight
						match self.graph.find_edge(prev, node) {
							Some(edge) => {
								*self.graph.edge_weight_mut(edge).unwrap() += 1;
							}
							None => {
								// if prev.index() != head.index() && prev.index() != node.index() {
								if prev.index() != node.index() {
                                    // 判断不为头节点和当前节点
									self.graph.add_edge(prev, node, 1);
								}
							}
						}
						prev = NodeIndex::new(*p);
					}
					node_indices.push(node); // 添加节点索引到向量
					i += 1;
				}
				AlignmentOperation::Ins(None) => {
					let node = self.graph.add_node(seq[i]);
					node_indices.push(node); // 添加节点索引到向量
					if edge_not_connected {
						self.graph.add_edge(prev, node, 1);
					}
					prev = node;
					edge_not_connected = true;
					i += 1;
				}
				AlignmentOperation::Ins(Some(_)) => {
					let node = self.graph.add_node(seq[i]);
					node_indices.push(node); // 添加节点索引到向量
					self.graph.add_edge(prev, node, 1);
					prev = node;
					i += 1;
				}
				AlignmentOperation::Del(_) => {} // we should only have to skip over deleted nodes and xclip
				AlignmentOperation::Xclip(_) => {}
				AlignmentOperation::Yclip(_, r) => {
					i = *r;
				}
			}
		}
		node_indices // 返回包含节点索引的向量
	}
	

    fn correct_base(&self, node: NodeIndex<usize>, base: u8) -> u8 {
        let total_weight: i32 = self.graph.edges(node)
            .map(|edge| *edge.weight())
            .sum();
        let node_weight = *self.graph.node_weight(node).unwrap();
        if total_weight > 0 && (node_weight as f32 / total_weight as f32) < 0.1 {
            node_weight
        } else {
            base
        }
    }
}

pub fn get_node_weights(graph: &POAGraph, node_indices: &Vec<NodeIndex<usize>>) -> Vec<u8> {
	node_indices.iter()
		.map(|&index| graph[index])
		.collect()
}

#[derive(Default, Clone, Debug)]
pub struct Aligner<F: MatchFunc> {
	pub first_seq_index: Vec<NodeIndex<usize>>,
	traceback: Traceback,
	query: Vec<u8>,
	poa: Poa<F>,
}
impl<F: MatchFunc> Aligner<F> {
    /// Create new instance.
    pub fn new(scoring: Scoring<F>, reference: TextSlice) -> Self {
		let (poa, first_sequence_indices) = Poa::from_string_get_index(scoring, reference);
        Aligner {
			first_seq_index: first_sequence_indices,
            traceback: Traceback::new(),
            query: reference.to_vec(),
            poa: poa,
        }
    }
	pub fn alignment(&self) -> Alignment {
        self.traceback.alignment()
    }
	pub fn add_to_graph(&mut self) -> &mut Self {
        let alignment = self.traceback.alignment();
        self.poa.add_alignment(&alignment, &self.query);
        self
    }
	pub fn add_to_graph_and_correct(&mut self) -> Vec<NodeIndex<usize>>{
		let alignment = self.traceback.alignment();
		let index: Vec<NodeIndex<usize>> = self.poa.add_alignment_get_index(&alignment, &self.query);
		index
	}
	pub fn global(&mut self, query: TextSlice) -> &mut Self {
        // Store the current clip penalties
        let clip_penalties = [
            self.poa.scoring.xclip_prefix,
            self.poa.scoring.xclip_suffix,
            self.poa.scoring.yclip_prefix,
            self.poa.scoring.yclip_suffix,
        ];

        // Temporarily Over-write the clip penalties
        self.poa.scoring.xclip_prefix = MIN_SCORE;
        self.poa.scoring.xclip_suffix = MIN_SCORE;
        self.poa.scoring.yclip_prefix = MIN_SCORE;
        self.poa.scoring.yclip_suffix = MIN_SCORE;

        self.query = query.to_vec();
        self.traceback = self.poa.custom(query);

        // Set the clip penalties to the original values
        self.poa.scoring.xclip_prefix = clip_penalties[0];
        self.poa.scoring.xclip_suffix = clip_penalties[1];
        self.poa.scoring.yclip_prefix = clip_penalties[2];
        self.poa.scoring.yclip_suffix = clip_penalties[3];

        self
    }
    pub fn graph(&self) -> &POAGraph {
        &self.poa.graph
    }
	pub fn consensus(&self) -> Vec<u8> {
        let mut consensus: Vec<u8> = vec![];
        let max_index = self.poa.graph.node_count();
        let mut weight_score_next_vec: Vec<(i32, i32, usize)> = vec![(0, 0, 0); max_index + 1];
        let mut topo = Topo::new(&self.poa.graph);
        // go through the nodes topologically
        while let Some(node) = topo.next(&self.poa.graph) {
            let mut best_weight_score_next: (i32, i32, usize) = (0, 0, usize::MAX);
            let neighbour_nodes = self.poa.graph.neighbors_directed(node, Incoming);
            // go through the incoming neighbour nodes
            for neighbour_node in neighbour_nodes {
                let neighbour_index = neighbour_node.index();
                let neighbour_score = weight_score_next_vec[neighbour_index].1;
                let edges = self.poa.graph.edges_connecting(neighbour_node, node);
                let weight = edges.map(|edge| edge.weight()).sum();
                let current_node_score = weight + neighbour_score;
                // save the neighbour node with the highest weight and score as best
                if (weight, current_node_score, neighbour_index) > best_weight_score_next {
                    best_weight_score_next = (weight, current_node_score, neighbour_index);
                }
            }
            weight_score_next_vec[node.index()] = best_weight_score_next;
        }
        // get the index of the max scored node (end of consensus)
        let mut pos = weight_score_next_vec
            .iter()
            .enumerate()
            .max_by_key(|(_, &value)| value.1)
            .map(|(idx, _)| idx)
            .unwrap();
        // go through weight_score_next_vec appending to the consensus
        while pos != usize::MAX {
            consensus.push(self.poa.graph.raw_nodes()[pos].weight);
            pos = weight_score_next_vec[pos].2;
        }
        consensus.reverse();
        consensus
    }
}

#[test]
fn poa(){
	let x = b"AGGTGGCTCAGCGCGGACAGCGCGTGGGGTAGGTCGTCCAGGCGCTCCACGGCGAGGCTCAGCGCGTCCGCCACCTTCTGGCCGTGGGCTCTGACTTGTGAGGAGCCGGGGCTCAGGTCCAGGTGGGAGAAGTAGGTCTTCGTGGCGGGGAAAGCCAGGAAGGTCCTGCGAGAGGAACGTGGGTGAGGCGCGTCCGGGGGGCTTGGGGAGGGAGGGCCCCTGGGGGCGGGGGCGCCCAGCCTGCCGCACCTTTCCAGGGCCTCTGTCGTGTAGACGCCGACGTTGCTGCCCAGCTTCTTCCACAGGGCGCGCACCAGCGCCCGGTCCTCCGCGGACAGCGCCATCCCCGCGCTGGAACCCGCCGCGCCCTGCGACCCCGCTGCGCCTGCGCGGACCGGCCAGGGGTCCCGCGGCCTCCAGGGGGCGCCGCGCCGCAGCTTCGCCGCCCGGGGGTCCAGGGCCTCGTCGGAGATCCCGGGACCGCGCGCGCTCAATCTCCGACAGCTCCGACCGGGGCCCTCTGAATCGCCCTGGATCCAAGGCCCCCCGTCTTAGGTGTTCACTATATGGTATTTAGTCTCTACTAAAAACACAAAAATTGGCCGGGTGCGGTGAGTGCAGTGGCGCCATTCCGGCTCACTGCAACCTCCGCCTCCCGGGTTCAAACGATTCTCCTGCCTCAGCCTCCCGAGTAACTGGGATCACAGGTGCATGCCACCATGCCAATTTTGTATTTTTAGTAGAGACAGGGTTTCACCGTGTTAGGCTGATCTCAAACTCCCGACCTCAGGTGATCCTCCCGCCTCAGCCTCCCAAAGTTCTGGGATTACAGGCGTGAACCACCGCGCCCGGCCATGTACATTATTATAGCACTGAATTTCGGCATGGATTTGGCATTGGAAGTCAGTTTATTCTGCCTCTTGATTCTATTGTGACACTTGCCCGAGTCTGTTTTCTCATTTATTATATAAAATGAGGGTGACAGTTTCTTCTCAGATTGTGAGAATTAAACGGATAACAACTGTGAGGTGCCTGACTCGTAGCAGGTGTTTCTTCAGGGCAGTGAACATTTATTGAGTACCTGGGATGTAACCCAAGAATAATAGTAAATGCATATGCAGCACTCACTCTGCTGTACCCCTTATATCTTTTTTGTTTTAAAAAAGGTGATAAAATAGACATAAATTACCATTGTAGCCATTTTTCAGGTATAAAATTTGGTGATATTAAGAACACTTAGTGTTATGCAGCCATCACCACTATCTAGGTCAGAACTTTTTCATCATCCCAAACGCCACCTCGGTACTCATTAACAGTAACTCTGGCCCCGGCAGCCACTAATCCGCTTTTGGTCCCTATGGATTTGGCTGTAGCGGATATTGCATATAGATGGGGTCATGCAGTTTGTGACCTTTGGTGTCCGGCTTCTTTTACTTGTCAAATATTTTCAGGGTGCAGCTACATCCTATCAGTGCTTCATTTCTTTTTGTGGCCAAATAATATTCCATTGTGTGTATAGACCACCACATTTTGTTTACCCATTCATCTACTGATGTATAATTTTTATAATTTTTTTTGAGATGGAGTTTTGCTCTTGTTGCCCAGGCTAGAGTGCAATGGTGTGATCTCGGGTCACTGCAACCTCTGCCTCCTGGGTTCAAGTGATTCTCCTACCTCAACCTCCCGAGTAGCTGGGATTACAGGCGTGCGCCACCACACCCGGCTAATTTTTTGTATTTTTAGTAGAGATGGGGTTTCCCCATATTGGCCAGGCTGGTCTCGAATTCCTGACCTTGTGATCTGCCCACCTTGGCTTCCCAAAGTGCTGGGATTACAGGCGTGAGCCACTGTGCCCAACCACCGTGCTAATTTTTAAATTAAATTAAATTAGATTAAATTTTTTTTTGAGACAGAGTTTCGCTCTTGTTGCCCAAGCTGGAGTGCGATGGCGCGATCTTGACTCACTACAACCTCCGCCTCCCAAAGTGCCAGGTTTACAGGCATGAGTCATCACACCTGGACAATTTTTGTCTTTTTAGTAGAGACAGGGTTTCACTATGTTGGCTAGGCTGGTCTCGAACCCCTGACTTCAGGTGATCCACCTGCCTAAGCCTGGCAGAGTGCTGCAATTACAGGTGTGAGCCACCGCACTTGGCCAATTTTTATAATTTTTCTTAGGCTTAGTCAGCTGAGCAGCAGGAGTGGAGAAAGGAAATGATAAATTTATAACTGGTTCTGATCAGTGTGGTGTAAAGCGCTGCCTCAACCAGCCCTCTGCTGTACATTTGACATATTTCTACCTTTAGGTCAGACCTCCATTCTTCCTGATGAAATAACGAAGGCACCGTCTTCCTGACAGTTGATACTGTACCAGGCTGCAGGAAGGACATCTGAAGAGCACGGAGGCAGAGGGCCTGCATAAGTGTTTACCAAAAGTGATTCATGGGCTCAGGGAATAAAAGATATTCCTATCAGTTGAGGGCCAGGGCTTTATGGCTACTGCTCCGGGAACCCTCCCTCTGATACCCCAGGCGTGGCCAGCCTCAGAGCTGAGTGGGTGCTGAGCTGGGCCCTGAGCTGGGCCCAGGAGACAACGCCTGATCTTGACAGCCCCAGGGCGTCAGTGTCTGGTGTGTGAGCATTCTGGGGTGACCTGGGGGTGAAAGTAGAGCCCCTATTGAATTGCCTGGGGCCTACAGTTGCACCCCACACCCTGATTTCATCTAGACTGTGCTGTATCTGGGAGGTTTACGGAGACTAGTGGGAACGATGGGGGATCAATGCAGCCGGCATGCTGGAGTGGGACTTCTCTGACCTGCACACACCACCCATATGTGGGAGGACCTCTGGGAAGCCAGTGACGTCCCAGGCAGAAAGCCAGCCAGTTCTTGGCACAAGAGTGCTAAGAAACGTAAAGATGGATAAATTCATTATTTCAGCAAACCTGCATTGAATCTGAAAAGTCTGGGAATAAAACTCAGGAAAGAGCAAATGCATCCTCAAAGCACTCTAGGGTCCAGCGTTTTTCCAGCTGGACTTCGCAGTGGCTCCACTTTCCCTCCTCCATCCCCTCCTCCCGCCCCTGCCTTTTCCTACCCTATCAGCTGCAGAGAGGTTCTAGCCATGTGTGTCCCAGCTGCTGTCCACGCCCATGCCTGGCCTTGGGGAAAAAGCTCAGGCACCCCCAGGCTGCCGCCCACTCAGACTTTATTCCGACCACGGGGGTACGGGTGCAGGAAGGGGAGGAGGGCTGGGGGGAGGAAGCAGGGGCAGAAGCATGGCCACCGGGCTCCCAGCAACAGTATTTGGAGGTCAGCACGGTGCTCACAGAAGCCCAGGAACTTGTCAGAGGCGTGCACCGCAGGAGTGAACTCGGCGGGGGAGGTGGGCGGCCAGGGTCACCAGCAGGCAGTGGCTTGGGAGCTGTGCAGAGAAGAGGGTCAGTGGGGCCGGGCCCAGGCCCGCAGCCGCCGCCTGCGCTTCACCTCCCGCAACCCGCGTGATCTGCCCTGCGAGGAAGGCGCCATCTCGCCCCTCGACGATCGCTCCCGGCCCCGCCGCTCACCTTGAAGTTGACCGGGTCCACCCGAAGCGCCTGTGCGCGTGCAGGTCGCTCAGGGCGGACAGCGCGTTGGGCATGTCGTCCACGTGCGCCACGGCGTTGGTCAGCGCGTCGGCCACCTTCTTGCCGTGGCCCTTAACCTGGGCAGAGCCGTGGCTCAGGTCGAAGTGCGGGGAGAAGTAGGTCTTGGTGGTGGGAAGGACAGGAACATCCCTGCGGGGAGAAGCAGGTGAGGGGTGGGGTTTTGGGTCGGGGGGGCCAGGACGGTTGAGGGGTGGCCTGTGGGTCCGGGCGGGCGAGAGCCCGGGTCCCGGGCAGGGAGGAGCCTCACCTCTCCAGGGCCTCCCCACCCATATAAAGCCAGCGTGCGCGCCCGACCTTACCCCAGGCGGCCTTGACGTTGGTCTTGTCGGCAGGAGACAGCACCATGGTGGGTTCTCTGAGTCTGTGGGGACCAGAAGAGTGCCGGGCCGCGAGCGCGCCAGGTTTATGCTTGGGGCGCGGGGGCACGCCCGGCCAGCCGGGCGGCACTCGTGGCTGGCGCGGAGCCCACCGGGGCGCGGCCTGGACCGCAGGGGAGTCCCGGGCGGGGCCTGCGCGGGGCCGGGCAGCCGCTCGGCCTCCCCTGGGGGTGCACGCGGGGCCGGGGGCCAGGACGTCTCACCCTCCACCCGCCACTCACTCCCGCCCATCCCGCTCGCCCGGGAACGGGGCGGTGGCGGCGCTGGCGGGGCCGGCCGTTGGGGTCGGGCGCTGTCGGCTCGTGCACCCCGGAGCGCAGCGCCGCCCTTTCCTTTCGGGCGCGGAGCGTCCTAGCGAGGAGAAAGTCAGCCCGCACCCCCGCCCCCGGCCTGGCACGCGCTGGACGCGCATCGACTCCAGCGGGATCGGGGAACACACGGGCGAGCGAGTGCGAGCCGGGAGCTTCCGCCATCCTGGGGCGGAGAGGAATGCGCACCCGGACGCCCTGGCCCATAGGAACTCAAAAGAATTTCTGCGCAGAGCCCCTCCTGCTCTCCAGCCTCTCGCTTCCCTGCAGCTCCTTTTCCCTCTGGCGATAGTCACTAGTGTGCTTGAGTGACAGGCACCGAAGGAATAAGCACCAGGACGCAAAAAGCACGGGGCTGGGCTTGGGGCGGCTCCAGATTCCAGACTCCTCTACCCGGGGTTGTGGGGTGGCGCGTGGGGTGGTGGGTGCTGGCTAAGCCCCAAAGTCATGGACTCACAGTTTGTGTGACTGACTCTCCCACAGCTAGAGGTCGTGGTTCACTGTGACACCCACAGCTCTCCAAAATGCTATCTGATGGAAGACTTGCTAGGTAAATACTGGTTAAACAGGGTAAACAAAGCAATAGCTGGAACCGGCTGGAGCATTCAACCTCCTCTGGGAGGTAGGCAGTCCTCTAACCATCACACTAAGTACACACAGAGGGGTGCAGACCTCATCACAGGCCCTGAGGAAGTCTGGAAGTGTGGACGAGGCATTCAAAGAGGGACCAGGGTTAGTCTGAGGGGCTTCCTGGAGGAGGTGAGACTTAAGGATGTATTAGGTGGAGGAGGTGGAGGAGGAGGAGGCTTGGCTCGGGGAGGGGGAAAGGACAGAGGAGACAAGGGGCCCTGTGGTTGGAAAGTGTGTGGTGGCAGGAGGCGGGATCGCAGTAAGTCATTTCCTGAGGGTCTGGCAGAGAGATCTTGCTTTGAGGAGTGCATCAGGTCAGCACCCTTCAGCCTGCTCCCTGCCCTGTCTCCCCAGGCCCCCACTCAGGCCCTGCCTGGGTGTCCAGGCAAAGCCAGGGTAAACTGCCTTGGGCAGAGAAGAAGTAGCTCCGACCAGCCCAGCAATGCTTCCCCCGTTGGATCTTCTCATTTCCCACTCCTGTCTGCCCTCTTCTGACTCTGCCCACAGCCTGAAGGAAAAAAAAAAGATAAGATGTGGGGGCTTTTTTTTTTTTTTTTTTTGAGGCGGAGTTTCGCTGTTGTTTTCCAGGCTGGAGTGCAATGGTGTGATCTTAGCTCACTGCAACCTCTGCCTCCTTGGTTAAAATGATTCTCCTGCCTCAGCCTCCGGAGTAGGTGGGTTTACAGGCATGGGCCACCATGCCTGGCTGATTTTTGTATTTTTAGTATAGATGGGGTTTCTACATGTTGATCAGGGTGGTCTCAAACTCCCGACCTCAGGTGATCCTCCCGCCTCGGCCTCCCAAAGTGCTGGGATTACAGGCGTGAGCCACCGTGCCTGGCGGAATTTTCGTTTATTTTTAAAGTGTTTCTGCTGAAATAACAATGCTCTGGGCTGTGAGGCGCAGGAAAAGTGAGTGGTGGAGAGGAGGGAACTCACAACAATTCAGAGGAGTCCATGTCCAAGGAAAGCCCGGTGACAGGGTCCTCTCCACCCTCCCCAAAGCTCAATCCAGTACTGTGAGTGTGGAGGGGAGCTTGCTGTGGGGCAGGGACTGGGTTTGGAAGTGCTGCTCATGTTCTAGGCTGGCTGAGAGCCAGATGGGCTTGTGCAACTTCCCAGCATGGATTCCAGGACTTCAGCACCGGGGACAGGGACGTTTCAGAAGTTGCTGGCCTTTGAAGGCAGAAAGCACCAGGAGGTCCCAGAGGCCCTGGGTGGGGCAGTGGCTGGGAGCTCAGGGCTGGGCACACAGTGAGTTCTGAATGCACCTTACAAACTGCAGGAGACGTGCAGGGTAAGGCCACAAGGTGCTCTGCTGAAGAGTCTCCCAGGTGTGGGCTGCCGGCTGTGAGCATCCCAGCCCCTTCCTCCTTCCCCACTGACTGGAGCTGCTCTCAGAGGGTCCTGCATGACACAGGAAGTGTGGTAGGGCCAGCACCCTCCTCAGCCCACCTACCCTCTGGCAGGAAGGGGTGGGAATGAGAGAAATGTTCTGGCACCTGCACCTGCACTGGGGACAGCCTATTTTGCTTGTTTTGTTTCGTTTTGTTTTGATGGAGAGCGTATGTTAACTGGGTTCAACAGGGGAGATTTCAAAGGGGTGGCGATGCTGGGACTCCCTGGTAAGACCCTGGATGCCTCTCTCCTCCCCAGCTTCCTGAGCCACTGCCTGCAGGCCTGGCACCTCTCCCAGGACAGGATGGCCAGCACCTTCTCTTGAGCAGGCCCACCTCAGCTTCCCCTCCCATACTCCTGTAGTTTCTCCTCCCCCAGCCCCATACAGCTGCAGAGGAGGTCCTTGGTCTGAGACAGGTAAACACCTCCATTGTTGGCACATTCCGGGGCAGAGAGAACCAGGCACACACAGGCTGCCGCCCACTCAGACTTTATTCAAAGACCAGGAAGGGCCAGTGCAAGGAGAGGAGGACCCGTTAGGGAGGCCCAGCGGGCAGGAGGAACGGCTACCAGGCTCAGCTTGACTGTTGGAGGGTCAGCACGGTGCTCACAGAAGCCAGGAACTTGTCCAGGGAGGCGTGCACCGCAGGTGAGCTCGGCGGGGAAGGTGGGCACGGCCAGGGTCACCAGCAGGCAGTGGCTTAGGAGCTGTGCAGAGGAAGGGGGTCAGTGCGGCCCAGGCCCACTGACAGCGCCGCCTGCGCTACACCTCCCGCAACCCGCGTGAGTCCTCTGCCCTGAGAGGAAGGCGCCATCCTCGCCCTCGACCCAGATCGCTCCCGGCCCGCCGCTCACCTTGAAGTTGACCGGGTCCACCCGAAGCTTGTGCGCGTGCAGGTCGCTCAGGGCGGACAGCGCGTTGGGCATGTCGTCCACGTGCGCCACGGCGTTGGTCAGCGCGTCGGCCACCTTCTTGCCGTGGCCCTTAACCTGGGCAGAGCCGTGGCTCAGGTCGAAGTGCGGGGAAGTGGGTCTTGGTGGTGAGGGGGAAGGACGGGAACATCCTGCGGGGAGAAGCAGAGTGAGGGGTGGGGGTTTGGGTCGGGGCCAGGACGGTTGGTGGCCTGTGGGTCCGGGCTGAACGAGAGCCCGGGTCGGGTACAGGGGAGGAGCCTCACCTCTCCAGGGCCTCCGCACCATACTCGCCAGCGTGCGCGCTGACGCCTTACCCCAGGCGGCCTTGACGTTGGTCTTGTCGGCGGGAGACAGCACCATGGTGGGTTCTCTCTGAGTCTGTGGGGACAGAAGAGTGCCGGGCCGCGAGCGCGCCGGGGTTTATGCTTGGGGCGCGGGGGCACGCCGGGCGGCGCTCGTTGGCTGGCGCGGAGCCCGGGGCGCGGCCTGGACCGCCGGGGGAGTCCGGGCGGGGCCTGCGCGGGGCCGGGCGGCGGGCTCGGCCTCCCCTGGGGGTGCACGCGGGGCGGGGGGCCAGGACGTCTCCACCCTCCACCCCGCCACTCCACTCCCGCCCATCGCTCGCCCGGGAGCAGGGCGGTAGCGGCGCGGGGGCCGGCCCCGTTGGGGTCGGGCGCTGTCAGCTCGTGCACCCCGGGCGCAGCGCCGCCCTTTCCTTTCGGGCGCCGGAGCGTCCTAGCGAGGGAGAAAGTCAGCCCGCACCCCCGCCCCGGCCTGGCACGCTGGACGCGCATCGACTCAGCGGGATCGGGGAACACACGGGCGAGCGAGTGCGAGCCGGGAGCTTCGCCCAATCCTGGGGCGGAGAGAATGCGCGCACCCGGACGCCCTGGCCCATAGGAACTCAAAGGTGCTGGAGCCCTCCTGCTCTCCAGCCTCTCCTTCCCTGCAGCTCCTTTCCCTCTGGCGATAGTCACTAGTGTGCCAGGTGACAGGCACACGGGAAGGAACAAACACCAGGACGCAAAAGCACGGGGCTGGGCTGGGGGCGGCTCGAATTCAGACTCCCTCTACCCGGGGTTGTGGGGTGGCGCGCGTGGGGTGGTGGGTACTGGCTAAGCCCCAAGTCATGGACTCACAGTTTGTGTGACTGACTCTCCCACTGGCCCCTAGAGGTCGTGGTTCCTGTGACACCCCAGCTCTCCAAATGCTATCTGATGGAAGACTTGCTGGGTAAATGCTGGTTGTAAACAGATAAACAAACTTGGCTCTGGGTAGGGAAAGGACAGGGGAGAGGAAGGCCCCTGTGAGTTGCAGAATGTAGGTGAGGGGATGTGGGGTGAGGAAGGGAAGGGGTGGACTTGGCGAGGGGAGGGTGGGGGTGCAGGAGTGCCAGTGGGTCACATCCCCTGCTCATGAGAAAGTGATGGTTTGGCCAGGCTCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCAAGGCGGGTGGATCAGCTGAGGTCGGGAGTTCGAGACCAGGCTGACCAACATGGAGAAACCCCGTCTCTACTAAAAATACAAAAATTAGCTGGGTGTGGTGGTGCACGCCTGTAATCCCAGCTACTCGGGAGGTGGAGGCTGGAGAATCCTTTGAACCAGGGAGTTGGAGGTTGCAGTGAGCCGAGATCGTGCCACAGCACTCTAGCCTGACGACACAGCGAGACTCTGTCTCAAAAAAATAAATACATAAGGCCGGGCGCGGTGGCTCACACCTATAATCCCAGCACATTGGGAGGCTGAGGCGGGCAGATCACGAGGTCAGGAGATCGAGACCATCCTCGCTAATGCAGTGAAACTCTGTCTCTATTAAAAATAGAAAAAATTAGCCAGGCGTAGTGGGGGGCGCCTGTAGTCCCAGCTACTCAGCTACTCGGGAGGCTGAGGAGAGGAGAATGGCGTGAACCTGCGAAGGTGGAGCTTGCAGTGAGCCGAGATTGTGCCACTGCACTCCAGCCTGGGTGACAGAGCGAGACTCCGTCTCAAAAAAAAAAAAAAAAAAGTAATAATAATAATAAAAATAAAAAAATAAAGTCATGATTTTTTTGGGGGATCATGATGGAAACATAGTAATAATCAGTGAGACTGTGGAATGGCACAGAGCTGAATGAACTGGCTGAAAGGGATGCAGAGGAGACACTTCACTGAGAATAGGAAGTTGTACAGGTATCTACAGTATGATGGTATTTTTGTCCAAAGAAGGCCTGTGGCAGTAGTTGTAGATGTAGCTGTGTTCCCTCAGCATGGGGATAGGGCCCATGGGCTGAGTTCAAACCTATACATTTTGGTCTCTTCTGGGGTGGGTGTGAGGAGACAGGAAGAGACACTCTCCTGCCAGAAGTAACACAAAAAGGGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGAAAGCATGGATAGATTTTAGTTTTATTTATTTATTTATTTATTTATTTTTTGAGACGGAGTTTTGCTCTTGTTGCCCAGGCTGGAGTGCAATGGCGTGATCTCAGCTCACCTCAACCTCCACCTCCCGGGTTCAAGCGATTCTCCTGCCTCAGCCTTCCTAGTAGCTAGGATTACAGGCATGCGCCACCACACCCAGCTAATTTTGTATTTTTAGTATAACTGGGGTTTCTCCATATTGGTCAGGCTGGTCTCAAACTCCCGACCTCAGGTGATCCTCTCGCCTCAGCCTCCCAAAGTACTGGGATTACAGGTGTGAGCCACCGTGCCTGTCGGATTTTAGTTATTTTTTAAAAGTGTTTCTGTTGAAATAACAATGCTCTGGGCTGTGAGGCACAGGAAGAGTGGGTGGTGGGAGAGGAGGGAACTCCCCAACAATTCAGAGAGGTCCATGTCCAGAAAAGCGGTGACAGGGTCCTCTCCACCCCTACCCAAAGCTCAATTCAGTACTGTGAGTGTGGAGGAGGGAGGCTTGCTGTGGGGCAGGGGCTGGGCTGGAAGTGCTGCTCGTGTTCCTAGGCTGGCTGGGGAGCCAGATGGGCTGTGCAACTTCCCAGCATGGATTCCAGGACTTCAGCACCGGGGACAGGGACGTTTCCAGAAGTTGCTGGCTTGGCTTTGAAGCAGAAAGCACCAGGAGGTCCCAGAGGCCCTGGGTGGGGCAATTAGCTGGGAGCTCAGAGCTGGGCACACAGTGAGTTCTGAATGCACCTTACAAACTGCAGGAGACGTGCAGGAGTAAGTGCAAAGTGCTCTGCTGAAGAGTCTCCCCAGGTGTGGGCTGCCGGCTGTGAGCATCACAGCCCCTTCCCTCCTTCCCCACTGACTGGGCTGCTCTCAGAGGGTCACCATGACACAGGAAGTGTGAGTAGGGCCGGGCACCCTCACTCCAGCCACCTACCCTCTGGCAGGGGGTGAGAATGAGAGAAATGTTCTGGCACCCGTTTTCTACTTGCATTGGGGACAACCTATTTTGCTAGTTTGTTTTGTTTCGTTTTGTTTTGATGGGAGGCGTATGTTAACTGGGTTCAACCAGGAGATTTCAAAGGAGGGGTGGCGATGCAGGACTCACAGGTAGAACCCTGGATGCCTCTCTCCTCCCCAGCTTCTGAGCCACTGCCTGCAGGCCTGGCACCTCTCAGGACAGGGGATGGTTCAGCACCTTCTCTTGAGCAGGCCCACCTCGGAGGTCCTTGGTCTGGAACAGGTAAACACCTCCACTGTTGGCACATTCCGGGACAGAGAACCCAGGCACACACAGGCTGCTGCCCTGCTCGGACTTCATTCAAAGATCAGGAAGTGCTGGGGCAGGGGGGGGGGGGAAGGAGGAGGACACTGGTGGGGCCCAGCAGGCAGGGAGGCATGACCACCAATGCTCTAGCTTAATGGTGTTTAGGTCAACACAGCACTCACAAAGGCGCAAATAATTTGTCCATGGAGGCACCACATGGGTAAATGGGAAGGGTGTGGCAGGCCTAGGTCACCAGCAGGCATTTACCCAGGAACTGTGGAGAGAGAGTCAGGGGATGCAGATACTACTGCCACCTGTGCTCCACCTCCCTGAACGCGTGGGAGGGTCTGTCCTATGAGGAAGGCCAGCTGCCCCCACACCCTAGATTCCTCCCAGCCCAAAGCTCTCCCTGATGTTGCCTGGGTCCACCCACAGCTCGTGGACATGCAGCTTCCTCAGCTCAGACACATCATTGGGCATGTCATCTAAGTGGACCACGGCATTGGTCAGCGTGTCGGCCACCTACTTGCATTGTTTCTTGATCTGGGCTGAGCCATAACTCAGGGCCAAGTGAGAGGAAGGTCAGGAACATGCTTCAGGGGGAAGCAAAGAGTGAGGGGCTGGGGGCTGGGTCCAGGGTAAGAATAGTTAGGGGTGGCCAGTGGAGCCCTGGAGTGGGCGTTGGACCCAAGGTGTGAGCAGGGAGAGGAGGGTTCTTGCCTCTCCAGGGCCTCCGTGGCATAGCCAGCAGTGTGGTTGCCAACTTTCTCCCAGACCGCCTTGGTGTTAGCCTTGTCCTCAGGAGACAGCACCACGGCAGCTTATTTCTGAGTCTGTGTAAGCTCAAAAGTGAAAGGTGGGGGAGAGTAGGAGCTCTTATGCCTGGGATGCAGGGTCATTGGCTGGCAAGTTGCCAGGAGGTGCCTGGCACAGGAAAGGGGAGGAGGAGTGAGTGAGAGTAGACAGGTTGCCCAAGGCTGCCCTGCGTCTCTTAGTCCTGCCCAAACAGAGACCAGCTCTCCACATGGTAGCCTTCCCGCTGGCATTGTCCTGAGCTTCCACACAGGGTCCTGAGCAGCTTATACCCATAAGGTATCCACATATGGCCCAGACCCAGGCAGCTTCCTGGATGGGGAGAAGAGCACAGTTTTTCCAAAGTGAAACATAAGCTGCTTTCACACCAAACTCTGCCTGGGACAAAAAGCTGTTAATTATCCATGTAAACAGTTGCGGGGGCTGGGTGTGGTGGCTCACGCCTGTCATCCCAGCACTTTGGGAGACCGAGGAGGGCAGATCACTTGAGGTCAGGAGTTTGAGACCAGCCTGATCAACATGGTGAAACTGTGTCTCTATTAAAAATACAAAAATTAGCCAGGAATGGCGGCCCATGCCTGTAATCCCAGCTACACGGGAGACAGAGGCAGGAGAATCGCTTGAACCCAGGAGGTTGAGGCTGCAGTGAGCCAAAACTGTGCCACTGCACTCCAGCCGGGGTGTCAGAGCAAGGGCCTATCTCAAAAACAAAAACAAAGCAAAAAAAAAAAAAACCCCCAAAAACAAAAAAACAGTTGGGGGATGTAGATAACGTGGGGCCAAGTTCTCCAAACTACCGGGCTGGGGTCCCAGGCAGTGCAGGGCAACAGGTCACTGGCACTGACTGCTGGGGGGAGTTTGGCCGACTCACAGTCAGGGCTCCACTTCCAGCAGCTCTCTTCAGACCTCGTTGCTGACTGCCAGGTCAGACCTGGGTGCCAGGCCTGGTCCAGTGCTCCCCACACAGACCCAGGATGGCCAACCCAGAGCCAACCCTATTCAGGGTTTGGAGACCTCCAAGGTCCCACTCTTTCCTCAGCCCCTATTCTTTGTCTGAAAAGCCTGGGGTGGAGCCCTAAGATTCCCTTCCCAGTCCTGGGGAGGGAGAAGAGGGAACCTGGGGAAGGGAGTGCCTTGGCCTGGGCACCTGCAGAGATTGAGACACGGCAGCTGGGGCTGGGTGGAGGGGTGGGGTGGCAACGACCCTGCAGGCACCCACCCTGTGTTATGATTTCAAAAATTATTTCAAAAATCTCGTGAGTTGGACGTTGTTTTTCTATTCGATAGGTAAAGGCCTGAGAGCTCTGGGACACGCCCAGGGTCGCACAGCTTGTCATTGAGCAGGCTGGAAGACTAAACCTGCTCTGTCTGCCAAATTGTGGCCACTGAGGACTTCCTCCCGACCTCATCCC";
	let y = b"AGGTGGCTCAGCGCGGACAGCGCGTGGGGTAGGTCGTCCAGGCGCTCCACGGCGAGGCTCAGCGCGTCCGCCACCTTCTGGCCGTGGGCTCTGACTTGTGAGGAGCCGGGGCTCAGGTCCAGGTGGGAGAAGTAGGTCTTCGTGGCGGGGAAAGCCAGGAAGGTCCTGCGAGAGGAACGTGGGTGAGGCGCGTCCGGGGGGCTTGGGGAGGGAGGGCCCCTGGGGGCGGGGGCGCCCAGCCTGCCGCACCTTTCCAGGGCCTCTGTCGTGTAGACGCCGACGTTGCTGCCCAGCTTCTTCCACAGGGCGCGCACCAGCGCCCGGTCCTCCGCGGACAGCGCCATCCCCGCGCTGGAACCCGCCGCGCCCTGCGACCCCGCTGCGCCTGCGCGGACCGGCCAGGGGTCCCGCGGCCTCCAGGGGGCGCCGCGCCGCAGCTTCGCCGCCCGGGGGTCCAGGGCCTCGTCGGAGATCCCGGGACCGCGCGCGCTCAATCTCCGACAGCTCCGACCGGGGCCCTCTGAATCGCCCTGGATCCAAGGCCCCCCGTCTTAGGTGTTCACTATATGGTATTTAGTCTCTACTAAAAACACAAAAATTGGCCGGGTGCGGTGAGTGCAGTGGCGCCATTCCGGCTCACTGCAACCTCCGCCTCCCGGGTTCAAACGATTCTCCTGCCTCAGCCTCCCGAGTAACTGGGATCACAGGTGCATGCCACCATGCCAATTTTGTATTTTTAGTAGAGACAGGGTTTCACCGTGTTAGGCTGATCTCAAACTCCCGACCTCAGGTGATCCTCCCGCCTCAGCCTCCCAAAGTTCTGGGATTACAGGCGTGAACCACCGCGCCCGGCCATGTACATTATTATAGCACTGAATTTCGGCATGGATTTGGCATTGGAAGTCAGTTTATTCTGCCTCTTGATTCTATTGTGACACTTGCCCGAGTCTGTTTTCTCATTTATTATATAAAATGAGGGTGACAGTTTCTTCTCAGATTGTGAGAATTAAACGGATAACAACTGTGAGGTGCCTGACTCGTAGCAGGTGTTTCTTCAGGGCAGTGAACATTTATTGAGTACCTGGGATGTAACCCAAGAATAATAGTAAATGCATATGCAGCACTCACTCTGCTGTACCCCTTATATCTTTTTTGTTTTAAAAAAGGTGATAAAATAGACATAAATTACCATTGTAGCCATTTTTCAGGTATAAAATTTGGTGATATTAAGAACACTTAGTGTTATGCAGCCATCACCACTATCTAGGTCAGAACTTTTTCATCATCCCAAACGCCACCTCGGTACTCATTAACAGTAACTCTGGCCCCGGCAGCCACTAATCCGCTTTTGGTCCCTATGGATTTGGCTGTAGCGGATATTGCATATAGATGGGGTCATGCAGTTTGTGACCTTTGGTGTCCGGCTTCTTTTACTTGTCAAATATTTTCAGGGTGCAGCTACATCCTATCAGTGCTTCATTTCTTTTTGTGGCCAAATAATATTCCATTGTGTGTATAGACCACCACATTTTGTTTACCCATTCATCTACTGATGTATAATTTTTATAATTTTTTTTGAGATGGAGTTTTGCTCTTGTTGCCCAGGCTAGAGTGCAATGGTGTGATCTCGGGTCACTGCAACCTCTGCCTCCTGGGTTCAAGTGATTCTCCTACCTCAACCTCCCGAGTAGCTGGGATTACAGGCGTGCGCCACCACACCCGGCTAATTTTTTGTATTTTTAGTAGAGATGGGGTTTCCCCATATTGGCCAGGCTGGTCTCGAATTCCTGACCTTGTGATCTGCCCACCTTGGCTTCCCAAAGTGCTGGGATTACAGGCGTGAGCCACTGTGCCCAACCACCGTGCTAATTTTTAAATTAAATTAAATTAGATTAAATTTTTTTTTGAGACAGAGTTTCGCTCTTGTTGCCCAAGCTGGAGTGCGATGGCGCGATCTTGACTCACTACAACCTCCGCCTCCCAAAGTGCCAGGTTTACAGGCATGAGTCATCACACCTGGACAATTTTTGTCTTTTTAGTAGAGACAGGGTTTCACTATGTTGGCTAGGCTGGTCTCGAACCCCTGACTTCAGGTGATCCACCTGCCTAAGCCCTGGCAAAGTGCTGCAATTACAGGTGTGAGCCACCGCGCCCGGCCAATTTTTATAATTTTTTCTTAGGCTTAGTCAACTGAAGCAGCAGGAGTATGAAAGGAAATGATAAATTTATAACTGGTTCTGATCAATTAATTGTAAAGCTGCACTCCAACCAGCCCTCTGCTGTACATTTGACATATTTCTACCTTTAGGTCAGACCTCCATTCTTCTGATGAAATAACGAAGACACCGTCTTCCTGACAGTTGATACTGTACCCACGAACTGCGGGAAGGACATCACAAACGCAGGCAGAGGGCCTGCCATAGGTGTTTGCAAGGGTGATTCATGGGCTCAGGGAATAGCAAAGATATTCCTATCAGTTGAGGGCCAGGGCTTTATGGCTGGCTCCGGGAACCCTCCCTCTGATACCCCAGGCACTGCAAACTATGCTCAACAGAGAGAAAACACCCAGGGAGAAGGGACCTTCTAGCCAGGCCTCTGGCCCCGCTCTCCTCTCCAGTCGCAATGGGACAGGCTCCCCAGCCTCAGGGAGGCCTGAGTGGGTGCTGAGCTGGGCCCTGAGCTGGGCCCAGAGACAAACGCCTGATCTTGACAGCCCCAGGGCGTCAGTGTCTGGTGTGTGAGCATTCTGGGGTGACCTGGGGGTGAAATTAGAGCCCCTATTGAATTGCCTGGGGCCTGCAGTTGCACCCCGCACCCCTGATTTCATCTAGACTGTGCCTGTATCTGGGAGGTTTTACGGAGACTAGTGGGAACGATGGGATCAATGCAGCCATGCTGGAGTGGGACTTCTCTGACTTACCCACCACCCATATGTGGGAGGACCTCACAGAGCCCAGTGACGTCCCAGGCAGAAAGCCAGCCAGTTCTTAACAAGAGTGCCCAGAAACGTAAAGATGGATAAATTCATTATTTCAGCAAACCTGCATTGAATCTGAAAAGTCTGGGGAATGTCAGGGAAGAACACCCAAATGCATCCTCAAGCACTCTAGGGTCCAGCGTTTTTCCAGCTGGACTCGCAGTAGCTCCGCTTCCCTCCTCCATCCCCTCCTCCCGCCCCTGCCTTTTCCTACCCCTATCAGCTGCAGAGGTTCTAGCCATGTGTGTCCCAGCTGCTGTCCACGCCCATGCCTGGCACGTTTGCTGAGGGGAAAACTCAGGCACACACAGGCTGCGCCCACTCAGACTTTATTCAAAGACCACGGGGGTACGGTGCAGGAAGGGGAGGAGGGGCTGGGGAAGGCCCAGAGGCAGAAGCATGGCCACCGAGGCTCCAGCAGCGGTATTTGGAGGTCAACTTTGAGTGCTCACAGAAGCCAGGAACTTGTCCAGGGAGGCGTGCACCGCAGGGAGTGAACTCGGCGGGGAGGTGGGCGGCCAGGGTCACCAGCAGGCAGTGGCTTAGGAGCTGTGCAGAGAAGAGGGTCAGTGGGGCCGAGGGCCCAGGGCCCGCAGCCTGCGCCCTCCGCAGCCCGCGTGATCCTGCCCTGCGAGGAGGCTGCATCTCGCCCCTCGACCCAGATCACTCCCCGGCCCGCCGCTCACCTTGAAGTTGACCGGGTCCACCCGAAGCTTGTGCGCGTGCAGGTCGCTCAGGGCGGACAGCGCGTTGGGCATGTCGTCCACGTGCGCCACGGCGTTGGTCAGCGCGTCGGCCACCTTCTTGCCGTGGCCCTTAACCTGGGCAGAGCCGTGGCTCAGGTCGAAGTGCGGGAAGTAGGTCTTGGTGGTGGGGAAGGACAGGAACATCCCTGCGGGGAGAAGCAAGTGAGGGTGGGGTTTGGGTCGGGGCCAGGACGGTCGAGGGAGTGGCCTGTGGGTCCGGGCGGGCGAGGCCCGGGCCTCACCTCTCCAGGGCCTCCGCACCATACTCGCCAGCGTGCGCGCCGACCTTACCCCAGGCGGCCTTGACGTCTTGTCGGCAGGGGCAGCGCCATGGTGGGTTCTCTCTGAGTCTGTGGGGACCAGAAGAGTGCCGGGCCGCAGGCGCGCCAGGGTTTATGCTGGGGGCGCGGGGGCACGCCCGGCCGGGCGGCAGCTCGTGGCTGGCGCGCGGAGCCGGGGCGCGGCCTGGACCACCAGGGAGTCCCGGGCGGGGCCTGCGCGGGGCCGGGCGGCGGGCTCAGCCTCCCCCTGGGTGCACGCGGGGCGGAGGCCAGGACGTCTCCGCCCCTCCACCCGCCACTCCCGCCAGCGCACCCATCCCCACTCGCCGGGAGCAGGGCGGTAGCGGCACTGGCGGGGGCCGGCCCGTTGGGGTCGGGCACGCTGTCGGCTCGTGCACCCCGGGCGCAGCGCCACCCTTTCCTTTCGGGCGCCGGAGCGTCCTAGCGAGGGGCCACCCCGCCCGGCTAGCGCGCTGGACGCGCATCGACTCCAGCGGGATCGGGGAACACACGGGCGGAGGCGAGTGCGGAGGCAATCCTGGGGCGGAGAGGAATGCGCACCCGGACGCCCTGGCCCATAGGAACTCAAAAGAATTTCTGCGCAGAGCCCCTCCTGCTCTCCAGCCTCGCTTCCTGCAGCTCCCTTTCCCTCTGGCGATAGTCACTAGTGTGCTTGAGTGACAGGCACCGGGAAGGAATAAACACCAGGACGCAAAAAGCACGGGGCTGGGCTGGGGGCGGCTCCCAGATTCAGACTCCTCTACCCGGGGTTGTGGGGTGGCGCGTGGGGGTGGTGGGTGCTGGCTAAGCCCCAAGTCATGGACTCACAGTTTGTGTGACTGACTCTCCCACTGGCCTAGAGGTCGTGGTTCACTGTGACACCCCCAGCTCTCCAAATGCTATCTGATGGAAGACTTGCTAGGTAAATACTGGTTAAACAGGTAAACAAAGCAATAGCTGGAACCGGCTGGAGCATTCAACCTCCTGGGGAGGGTAGGCAGTCCTCTAACCATCACACAAGTACACACAGAGGTGCAGACCTCATCACAGGCCCAGTCTGGAGGTGTGGACAGGAGGCATTCAAGGGGACCAGGGTTAGTCTGAGGGCTTCCTGGAGGAGGTGAGACTTAAGTATGTATTAGGTGGAGGAGGTGGAGGAGGGAGGAGGCTTGGCTCGGGGAGGGGAAAGGACAGAGGAGACAAGGGCCCTGTGGTTGGAGAATGGAGGTGGGGAGAATGTTGGGGTGATGAGGCTGGCCAGCTCCTGCAGGGTGAGGGAAGGGAAGGGGTGGACTTCTGGGGAAGGGTGGGAGGTGCAGGAGTGCCAGTGGGTGCATCCCTGTGCATGAAAAGTGATGTTTTGGGGGGATGGTACTGAGGAGAAAACAGCCTGAGAAATCACTGATAAGTCATTTCCTGGGGGTCTGGCAGAAGATCTTGCTTTGAGTGCATCAGGTCAGCACCCTTCAGCCTGCTCCCTGCCCTGTCTCCCCCAGGCCCTTACTCGGGCCCTGCCCTGGGTGTCCAGGAGCAAGCCAGGGTGGCTGCCTTGGGCAGGGAAGTAACCGACCAGCTTAGCAATGCTTCCCCCGTTGGATCTTCTCATTTCCCCTCCCTGTCTGCCACCCTCTTCTGACTCTGCCCACAGCCTGAAGGAAAAAAAAGATAAGATGTGGGGGCTTTTTTTTTTGAGGCGGAGTTTCGCTGTTGTTTCCCAGGCTGGAGTGCAATGGTGTGATCTTAGCTCACTGCAACCTCTGCCTCCTTGGTTAAAATGATTCTCCTGCCTCAGCCTCCGGAGTAGGTGGGTTTACAGGCATGGGCCACCATGCCCGGCTGATTTTGTATTTTTAGTATAGATGGGGTTTCTACATGTTGATCAGGGTGGTCTCAAACTCCCGACCTCAGGTGATCCTCCCACCTCGGCCTCCCAAAGTGCTGGGATTACAGGCGTGAGCCACCGTGCCTGTCGGATTTTAGTTTATTTTTAAAGTGTTTTCTGCTGAAATAACAATGCTCTGGGCTGTGAGGCGCAGGAAGAGCGGGTGGTGGAGGGGAGAGAACTCCCCAACAATTCAGAGGGTCCATGTCCAGAAGGAAATTTGTCCTCTCCTCCCCAAAGCTCAATCCAGTACTGTGAGTGTGGAGGGGAGGCTTGCTGTGGGGCAGGGACTGGGTTGGAAGTGCTGCTCATGTTCCTAGGCTGGCTGAGAGCCAGATGGGCTGTGCAGCTTCCCAGCATGGATTCCAGGACTTCAGCACCACGGGGACAGGGACGTTTCCAGAAGTTGCTGGCAGCACAGGAAGCAGAAGCACCAGGAGGTCCCAGAGGCCCTGGGTGGGGCAATTAGCTGGGAGCTCAGGGCTGGGCACACAGTGAGTTCTGAATGCACCTTACAAACTGCAGGAGACGTGCAGGAGTGGGCCGCAGATGCTCTGCTGAAGAGTCTCCCCAGGTGTGGGCTGCCGGCTGTGAGGCATCCCAGCCCCTTCCTCCTTCCCCACTGACTGGGCTGCTCTCAGAGGGGTCACCGTGACTGGGAAGTGTGAGTAGGGCCAGCACCCTCACTCCAGCCACCTACCCTCTGGCAGGAGGAAGGGGTGGGAATGAGAGAAATGTTCTGGCATGCACTTGCTGGGGACAGCCTATTTGCTAGTTTGTTTTGTTTCGTTTTGTTTTGATGGAGAGCGTATGTTAACTGGGTTCAACCAGGAGATTTCAAAGGAGGGTGGCGATGCTGGGACTCCCTGGGTAGAACCTGGATGCCTCTCCTCCCCAGGCTTCCCTGAGCCACTGCCTGCAGAGCCCCTGGCACCTCTCAGGACAGGGGATGGTTCAGCACCTTCTCTTGAGCAGGCCCACCTCAGCTTCCCCTCCATACTCCCTGCAGTTCTCCCTCCCCAGCCCCATGCAGCTGCAGGGAGGTCCTTGGTCTGAGACAGGTAAACACCTCCATTGTTGGCACATTCCGGGACAGAGAGAACCCAGGCACCTGCCGCCCCACTCAGACTTTATTCAAGACCAGGAAGGGCCGGTGCAAGGGAGGGGAGGAGGGGCCCAGGGTGGGAGGCCCAGCGGGCAGGAAGGTGACTGCCAGGGCTCAGCTTGACGGTATTTGGAGGGCACGGTGCTCACAGAAGCCAGAACTTGTCCAGGGAGGCGGTGCACCTGGCGGGTGAACTCGGCGGGGAGGTGGGCGGCCAGGGTCACCAGCAGGCAGTGGCTTAGGAGCTGTGCAGAGGGTCAGTGCGGCCCAGGCCCGCAGCCACCGCCTGCGCTTCTCGCCCCTCGACCCAGATCGCTCCCGGCCCGCCGCTCACCTTGAAGTTGACCGGGTCCACCCGAAGCTTGTGCGCGTGCAGGTCGCTCAGGGCGGACAGCGCGCGTTGGGCATGTCGTCCGTGCGCCACGGCGTTGGTCAGCGCGTCGGCCACCTTCTTGCGTGGCCCTTAACCTGGGCAGAGCCGTGCTCAGGTCGAAGTGCGGGAGTAGGTCTTGGTGGTGGGAAGGACAGGAACATCCTGCGGGGAGAAGCAGAGTGAGGGGTGGGGTTTGGGTCCGGGGCCGGGACGGTTGGTGCCTGTGGGTCCGGGCGGGCGAGGAGCCCGGGTCCGGAGCGGGAGAGGCCTTACCTCTCAGGGCCTCCGCGCGCGGTGGGTCTGTGGGGACCAGAAGGTGCCGGGCCGCAGGCGCGCCGGGGTTTATGCTTGGGGCGCGGGGCTATGTGGGGCGCGGCCTGGCCGCGGGGAGTCCCGGAGCGGGGCCTGCGCGGGGCCGGGCGGCGGGCTCGGCCTCCCCTGGGGGTGCACGCGGGGCGGGGGCCCGGGACGTCTCCACCCTCCACCCGCCGCCTCCCCGCCCATCCCGCTCGCCCGGGAGGGGGGGGCGGCAGCTGGCGGGGCCGGCCCGTTGGGGTCGGGCGCTGTCGCTCGTGCACCCCGGGCCTGGCGCCGCCCTTTCCTTTCGGGCGCCGGAGCGTCCCTAGCGAGGGAGAAAGTCAGCCCGCACCCCCGCCCCGGCCCTGGCACGCGCTGGACGCGCCATCGACTCCAGCGGGATCGGGGGGAACACACGAGCGAGTGCGAGCCGGGAGGCCTCGCCCAATCCTGGGGCGGAGAGGAATGCGCCCACCCGGACGCCCTGGCCCATAGGAACTCAAAGAATTTCTGCGCAGAGCCCCCTCCTGCTCTCCAGCCTCGCTTCCTGCAGCTCCTTTCCCTCTGGCGATAGTCACTAGTGTGCTTGAGTGACAGGCACCGGGAAGGAACAAACACCAGGACGCAAAAAGCACGGGGCTGGGCTGGGGGCGGCTCCAGATTCAGACTCCTCTACCCGGGGTTGTGGGGTGGCGCGTGGGGTGGTGGGTGCTGGCTAAGCCCCAAGTCATGGACTCACAGTTTGTGTGACTCTCCCACTGGCCTAGAGGTCGTGGTTCACTGTGACACCCCCAGCTCTCCAAATGCTATCTGATGGAAGACTTGCTAGGTAAATACTGGTTGTAAACAGATAAACAAACTTGAGCTCCTGGGTAGGGAAAGGACAGGGAAGGCCCCTGTGGTTGCAGAATGTAGGTGAGAGGGGATGTGGGGTGAGGAAGGAAGGGGTGGACTTGGCGAGGGGAGGGTGGGAGGTGCAGGAGTGCCAGTGGGTGCATCCCTGCTCATGAAAAGTGATGGTTTGGCCGGGCTCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCAAGGCGGGTGGATCAGCTGAGGTCGGGAGTTCGAGACCAGGCTGACCAACATGGAGAAACCCCGTCTCTACTAAAAATACAAAAATTAGCTGGGTGTGGTGGTGCACGCCTGTAATCCCAGCTACTCGGGAGGTGGAGGCTGGAGAATCCTTTGAACCAGGGAGTTGGAGGTTGCAGTGAGCCGAGATCGTGCCACAGCACTCTAGCCTGACGACACAGCGAGACTCTGTCTCAAAAAATAAATACATAAGGCCGGGCGCGGTGGCTCACGCCTGTAATCCCAGCACATTGGGAGGCTGAGGCGGGCAGATCACGAGGTCAGGAGATCGAGACCATCCTCGCTAATGCAGTGAAACTCTGTCTCTATTAAAAATAGAAAAAATTAGCCAGGCGTAGTGGGGGGCGCCTGTAGTCCCAGCTTCTCAGCTACTCGGGAGGCTGAGGGAGGAGAATGGCGTGAACCTGCGAGGTGGAGCTTGCAGTGAGCCGAGATTGTGCCACTGCACTCCAGCCTGGGTGACAGAGCGAGACTCCGTCTCAAAAAAAAAAAAAAAAAAATAATAATAATAATAAAAATAAAAAAATAAAAGTCATGATTTGCAGGATCATGATGGAAACACATAGTAATAATCAGTGAGACTGTGGAATGGCACAGAGCTGAATGAACTGGCTGAAAGGGATGCAGAGGAGACACTTCACTGAGAATAGGAAGTTGTACACAGGTATCTACAGTATGATGGTATTTTGTCCAAAAAGAGAGCCTGTGGCAGTAGTTGTAGATGTAGCTGTGTTCCCTCAGCATGGGATGGGAGGCCCTGGGCTTGAGTTCCCAAACCCTTCATTTTGGTCTCTTCTGGGGGTGTGAGGGGCAGAAGAGACACTCTCCTACTTTAAGTAACACAAAAAAGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGAAAGCATGGATAGATTTTTAGTTTATTTATTTATTTATTTATTTATTTTTTGAGACGGAGTTTTGCTCTTGTTGCCCAGGCTGGAGTGCTATCTCAGCTCACCTCAACCTCCACCTCCCGGGTTCAAGCGATTCTCCTGCCTCAGCCTTCCTAATAGCTAGGATTACAGGCATGCGCCACCACACCCAGCTAATTTTGTATTTTTAGTATAACTGGGGTTTCTCCATATTGGTCAGGCTGGTCTCAAACTCCCGACCTCAGGTGATCCTCTCGCCTCAGCCTCCCAAAGTACTGGGATTACAGGTGTGAGCCACCGTGCCTGTCCGGATTTTAGTTTATTTTTTAAAGTGTTTCTGTTGAAATAACAATGCTCTGGGCTGTGAGGCACAGAGAGAGTGGGTGGTGGAGAGGGGAGGGAACTCCCCAACAATTCAGGAAGGTCCATGTCCAGAAGAAAAGCGGTGACAGGGTCCTCTCCACCCCTACCCAAAGCTCAATTCAGTACTGTGAGTGTGGAGGGGAGGCTTGCTGTGGCAGGGGCTGGGCTGGAAGTGCTGCTCGTGTTCCTAGGCTGGCTGAGAGCCAGATGGGCGTGCAACTTCCCAGCATGGATTCCAGGACTTCAGCACCGGGGACAGGGACGTTTCCAGAAGTTGCTGGCAGCACAGGAAGCAGAAGCACCAGGAAGTCCCAGAGGCCCTGGGTGGGGCAATTAGCTGGGAGCTCAGGGCTGGGCACACAGTGAGTTCTGAATGCACCTTACAAACTGCAGGAGACGTGCAGGAGTAAGGCCGCAAGGTGCTCTGCTGAAGAGTCTCCCCAGGTGTGGGCTGCCGGCTGTAGGCATCCAGCCCCTTCCTCCTTCCCCACTGACTGGGCTGCTCTCAGAGGGTCACCATGACACAGGAAGTGTGAGTAGGGCCAGCACCCTCACTCCAGCCCTCTGGCAGGAAGGGGTGAGAATGAGAGAAATGTTCTGGCACCTGCACTTGCATTGGGGACAGCCTATTTTGCTAGTTTGTTTTGTTTCGTTTTGTTTTGATGGAGAGCTGGGTTCAACCAGGGAGATTTCAAAGGAGGGTGGCGATGCTGGGACTCACTGGGTAAGACCCCTGGATGCCTCTCTCTCCTCCCCAGCCTACACAGGCCACTGCCTGCAGGCCTGGCACCTCTCAGGACAGGGGATGGTTCAGCACCTTCTCTTGAGCAGGCCCACCTCAGCTTCCCCTCCCATACTCCCTGCAGTTCTCCCTCCCCAGCCCCATGCAGCAGGAGGTCCTTGGTCTGAGACAGGTAAACACCTCCTGTTGGCACATTCCGGGACCCAGAGAACCCAGGCACACACAGGCTGCTGCCTACTCGGACTTCATTCAAAGATCAGGAAGTGCTGGGGCAGGGGAGGAGAACAGACAGGAGGCCCAGCAGGCAGGAGCATGACCACCAATGCCTCTAGCTTAATGGTGTTTTTAAGGTCCAACACAGCACTCACAAAAGCAAATAATTTGTCCATGGAGGCACCACATGGGTAAATGGGAAGGGTGTGGCAGGCCTAGGTCACCAGCAGGCATTTACCCAGGAACTGTGGAGAGAGAGTCAGGGGATGCAGATACTACTGCCACCTGTGCTCCACCTCCCTGAACGCGTGGGAGGGTCTGTCCTATGAGGAAGGCCAGCTGCCCCCACACCCTAGATTCCTCCCAGCCCAAAGCTCTCCCTGATGTTGCCTGGGTCCACCCACAGCTCGTGGACATGCAGCTTCCTCAGCTCAGACACATCATTGGGCATGTCATCTAAGTGGACCACGGCATTGGTCAGCGTGTCGGCCACCTACTTGCATTGTTTCTTGATCTGGGCTGAGCCATAACTCAGGGCCAAGTGAGAGGAAGGTCAGGAACATGCTTCAGGGGGAAGCAAAGAGTGAGGGGCTGGGGGCTGGGTCCAGGGTAAGAATAGTTAGGGGTGGCCAGTGGAGCCCTGGAGTGGGCGTTGGACCCAAGGTGTGAGCAGGGAGAGGAGGGTTCTTGCCTCTCCAGGGCCTCCGTGGCATAGCCAGCAGTGTGGTTGCCAACTTTCTCCCAGACCGCCTTGGTGTTAGCCTTGTCCTCAGGAGACAGCACCACGGCAGCTTATTTCTGAGTCTGTGTAAGCTCAAAAGTGAAAGGTGGGGGAGAGTAGGAGCTCTTATGCCTGGGATGCAGGGTCATTGGCTGGCAAGTTGCCAGGAGGTGCCTGGCACAGGAAAGGGGAGGAGGAGTGAGTGAGAGTAGACAGGTTGCCCAAGGCTGCCCTGCGTCTCTTAGTCCTGCCCAAACAGAGACCAGCTCTCCACATGGTAGCCTTCCCGCTGGCATTGTCCTGAGCTTCCACACAGGGTCCTGAGCAGCTTATACCCATAAGGTATCCACATATGGCCCAGACCCAGGCAGCTTCCTGGATGGGGAGAAGAGCACAGTTTTTCCAAAGTGAAACATAAGCTGCTTTCACACCAAACTCTGCCTGGGACAAAAAGCTGTTAATTATCCATGTAAACAGTTGCGGGGGCTGGGTGTGGTGGCTCACGCCTGTCATCCCAGCACTTTGGGAGACCGAGGAGGGCAGATCACTTGAGGTCAGGAGTTTGAGACCAGCCTGATCAACATGGTGAAACTGTGTCTCTATTAAAAATACAAAAATTAGCCAGGAATGGCGGCCCATGCCTGTAATCCCAGCTACACGGGAGACAGAGGCAGGAGAATCGCTTGAACCCAGGAGGTTGAGGCTGCAGTGAGCCAAAACTGTGCCACTGCACTCCAGCCGGGGTGTCAGAGCAAGGGCCTATCTCAAAAACAAAAACAAAGCAAAAAAAAAAAAAACCCCCAAAAACAAAAAAACAGTTGGGGGATGTAGATAACGTGGGGCCAAGTTCTCCAAACTACCGGGCTGGGGTCCCAGGCAGTGCAGGGCAACAGGTCACTGGCACTGACTGCTGGGGGGAGTTTGGCCGACTCACAGTCAGGGCTCCACTTCCAGCAGCTCTCTTCAGACCTCGTTGCTGACTGCCAGGTCAGACCTGGGTGCCAGGCCTGGTCCAGTGCTCCCCACACAGACCCAGGATGGCCAACCCAGAGCCAACCCTATTCAGGGTTTGGAGACCTCCAAGGTCCCACTCTTTCCTCAGCCCCTATTCTTTGTCTGAAAAGCCTGGGGTGGAGCCCTAAGATTCCCTTCCCAGTCCTGGGGAGGGAGAAGAGGGAACCTGGGGAAGGGAGTGCCTTGGCCTGGGCACCTGCAGAGATTGAGACACGGCAGCTGGGGCTGGGTGGAGGGGTGGGGTGGCAACGACCCTGCAGGCACCCACCCTGTGTTATGATTTCAAAAATTATTTCAAAAATCTCGTGAGTTGGACGTTGTTTTTCTATTCGATAGGTAAAGGCCTGAGAGCTCTGGGACACGCCCAGGGTCGCACAGCTTGTCATTGAGCAGGCTGGAAGACTAAACCTGCTCTGTCTGCCAAATTGTGGCCACTGAGGACTTCCTCCCGACCTCATCCC";
	let p = b"AGGTGGCTCAGCGCGGACAGCGCGTGGGGTAGGTCGTCCAGGCGCTCCACGGCGAGGCTCAGCGCGTCCGCCACCTTCTGGCCGTGGGCTCTGACTTGTGAGGAGCCGGGGCTCAGGTCCAGGTGGGAGAAGTAGGTCTTCGTGGCGGGGAAAGCCAGGAAGGTCCTGCGAGAGGAACGTGGGTGAGGCGCGTCCGGGGGGCTTGGGGAGGGAGGGCCCCTGGGGGCGGGGGCGCCCAGCCTGCCGCACCTTTCCAGGGCCTCTGTCGTGTAGACGCCGACGTTGCTGCCCAGCTTCTTCCACAGGGCGCGCACCAGCGCCCGGTCCTCCGCGGACAGCGCCATCCCCGCGCTGGAACCCGCCGCGCCCTGCGACCCCGCTGCGCCTGCGCGGACCGGCCAGGGGTCCCGCGGCCTCCAGGGGGCGCCGCGCCGCAGCTTCGCCGCCCGGGGGTCCAGGGCCTCGTCGGAGATCCCGGGACCGCGCGCGCTCAATCTCCGACAGCTCCGACCGGGGCCCTCTGAATCGCCCTGGATCCAAGGCCCCCCGTCTTAGGTGTTCACTATATGGTATTTAGTCTCTACTAAAAACACAAAAATTGGCCGGGTGCGGTGAGTGCAGTGGCGCCATTCCGGCTCACTGCAACCTCCGCCTCCCGGGTTCAAACGATTCTCCTGCCTCAGCCTCCCGAGTAACTGGGATCACAGGTGCATGCCACCATGCCAATTTTGTATTTTTAGTAGAGACAGGGTTTCACCGTGTTAGGCTGATCTCAAACTCCCGACCTCAGGTGATCCTCCCGCCTCAGCCTCCCAAAGTTCTGGGATTACAGGCGTGAACCACCGCGCCCGGCCATGTACATTATTATAGCACTGAATTTCGGCATGGATTTGGCATTGGAAGTCAGTTTATTCTGCCTCTTGATTCTATTGTGACACTTGCCCGAGTCTGTTTTCTCATTTATTATATAAAATGAGGGTGACAGTTTCTTCTCAGATTGTGAGAATTAAACGGATAACAACTGTGAGGTGCCTGACTCGTAGCAGGTGTTTCTTCAGGGCAGTGAACATTTATTGAGTACCTGGGATGTAACCCAAGAATAATAGTAAATGCATATGCAGCACTCACTCTGCTGTACCCCTTATATCTTTTTTGTTTTAAAAAAGGTGATAAAATAGACATAAATTACCATTGTAGCCATTTTTCAGGTATAAAATTTGGTGATATTAAGAACACTTAGTGTTATGCAGCCATCACCACTATCTAGGTCAGAACTTTTTCATCATCCCAAACGCCACCTCGGTACTCATTAACAGTAACTCTGGCCCCGGCAGCCACTAATCCGCTTTTGGTCCCTATGGATTTGGCTGTAGCGGATATTGCATATAGATGGGGTCATGCAGTTTGTGACCTTTGGTGTCCGGCTTCTTTTACTTGTCAAATATTTTCAGGGTGCAGCTACATCCTATCAGTGCTTCATTTCTTTTTGTGGCCAAATAATATTCCATTGTGTGTATAGACCACCACATTTTGTTTACCCATTCATCTACTGATGTATAATTTTTATAATTTTTTTTGAGATGGAGTTTTGCTCTTGTTGCCCAGGCTAGAGTGCAATGGTGTGATCTCGGGTCACTGCAACCTCTGCCTCCTGGGTTCAAGTGATTCTCCTACCTCAACCTCCCGAGTAGCTGGGATTACAGGCGTGCGCCACCACACCCGGCTAATTTTTTGTATTTTTAGTAGAGATGGGGTTTCCCCATATTGGCCAGGCTGGTCTCGAATTCCTGACCTTGTGATCTGCCCACCTTGGCTTCCCAAAGTGCTGGGATTACAGGCGTGAGCCACTGTGCCCAACCACCGTGCTAATTTTTAAATTAAATTAAATTAGATTAAATTTTTTTTTGAGACAGAGTTTCGCTCTTGTTGCCCAAGCTGGAGTGCGATGGCGCGATCTTGACTCACTACAACCTCCGCCTCCCAAAGTGCCAGGTTTACAGGCATGAGTCATCACACCTGGACAATTTTTGTCTTTTTAGTAGAGACAGGGTTTCACTATGTTGGCTAGGCTGGTCTCGAACCCCTGACTTCAGGTGATCCACCTGCCTAAGCCTACCAAAGTGCTGCAATTACAGGTGTGAGCCACCGCGCCCGGCCAATTTTTATAATTTTTTCTTAGGCTTAGTCAACTGAAGCAGCAGGAGTGAAGAGCGTAAATAAATTTATAACTGGTTCTGATCAATTAATTGTAAAAAGCGCTGCACTCCAACCAGCCCTCTGCTGTACATTTGACATATTTCTACCTTTAGGTCAGACCTCCATTCTTCCTGATGAAATAACGAAGACACCGTCTTCCTGACAGTTGATACTGTACCCACAGGCTGCGGAAGGACATCACAAACGCAGGCCAAGAGGGGCCTGCCATAGGTGTTTACCAAGGGTGATTCATGGGCTCAGGGAATGAAAAGATATTCCTATCAGTTGAGGGCAGGGCTTTATGGCTTCCAGCTCCGAAGTCCCTCTGATACCCCAGGCACTGCAAGCCTGCCCAGCAGAGAAAACACACACCAGGGAGGGACCTTCTAGCCAAGGCCTCTGGCCCCGCTCTCACTCTCCAGTCGCAATGGGACAGGCTCCCCAGCCTCAGGGGAGCTGAGTGGAGTGCTGAGCTGGGCCCTGAAGCTGGGCCCAGAGACAAACGCCTGATCTTGACAGCCCCAGGGCGTCAGTGTCTGGTGTGTAGGCATTCTGGGGTGACCTGGGGGTGATAGGGCCCTATTGGATTGCCTGGGGCCTGCAGTTGCACCCCGCACCCCTGATTTCATCTAGACTGTGCCTGTATCTGGGAGGTTTTACGGAGACTAGTGGGAACGATAGGGGATCAATGCAGCCATGCTGGGCTTCTCTGACCCACCACCCATATGTGGGAGGACCTCTGAAACCAGTGACATCCCAGGCAAGCCGAGTGCCAAAACGTAAAGATGGATAAATTCATTATTTCAGCAAACCTACATTGAATCTGAAAAGTCTGGGAATAAAACTCGGGAAAGAGCAAATGCATCCTCAAAGCACTCTAGGGTCAGCGTTTTCCAGCTGGACTTCGCAGTGGCTCCACTTTCCCTCCTCCATCCCCTCCTCCCGCCCCTGCCTTTTCCTACCCTATCCAGCTGCAGACGGGTTCTAGCCATGTGTGTCCCCAGCTGCTGTCCGCCCATGCCTGGCACGTTTGCTGAGGAAAAAACTCAGGCACACTGGCTGCCGCCCACTCAGACTTTATTCAAAAAGACCACGGGGGTACGGGTGCAGGAAGGGGAGGAGGGGCTGGGGGGAGGCCCAAGGGGGCAAGAAGCATGGCCACCGAGGCTCCAGCTTAACGGTATTTGGAGGTCAGCACAGTTATGCAGAAGCCAGGAACTTGTCCAGGAGGCGTGCACACGGGGTGAACTCGGCAGGGAGGTGGGCGGCCGGTCACCCAGCAGGCAGTGGCTTAGGAGCGCAGAAGAGGGGTCAGTGGGGCCCGAGGGCCCAGAGCCCGCGGCCGCCGCCTGCGCTACGCCTCCGCGGCCGCGTGATCCTCTGCCCTGCGAGGAAGGCGCCATCTCGCCCCTCGACCCAGATCGCTCCCGGCCCGCCGCTCACCTTGAAGTTGACCAGGTCCACCCGAAGCTTGTGCGCGTGCAGGTCGCTCAGGGCGGACAGCGCGTTGGGCATGTCGTCCACGTGCGCCACGGCGTTGGTCCAGCGCGTCGGCCACCTTCTTGCCCGTGGCCCTTAACCTGGAGAGCCGTGGCTCAGGTCGAGTGCGGGAAGTAGGTCTTGGTGGTGGGGAAGGACAGGAACATCCTGCGGGAAGAGCAGAGTGAGGGGTGGGGTTGGGTCCGGGGCCAGGACGGTTGAGGGTGGCCTGTGGGTCCGGGCGGGCGAGGAGCCCGAGTCGGAGCAGGGGAGGGGCCTCACCTCTCCAGGGCCTCCGCACCATACTCGCCAGCGTGCGCCGACCTTACCCAGGCGGCCTTGACGTTGGTCTTGTCAGCAGACAGCACCATGGTGGGTTCTCTCTGAGTCTGTGGGGACCAGAAGAGTGCGGGCCGCGAGCGCGCCAGGGTTTGTACTTGGGGCGCGGGGGCGCCCGGCCGGGCGGCGCTCATTGGCTGGCGCGGAGCCCGGGGCGCGGCCTGGACCGCAGGGGGTCCGGGCGGGGCCTGCGCGGGGCCGGGCGGCAGGGCTCGGCCTCCCCTGGGGGTGCACGCGGGGCGGGGGCCAGGACGTCTCCACCCTCCACCCGCCACTCACTCCCGCCCATCCCTCGCCCGGGGAGCAGGGCCGTGGCGGCGCTGGCGGGGCCAGCAGGGGTCGGGGCGCTGTCGGCTCGTGCCCCGGAGCGCAGCGCCCCGCCCTTTCCTTTCGGGCGCCCGGAGCGTCCTAGCGAGGGAGAAAGTCGGCCCGCACCCCCCCGCCCCCGGCCTGGCACGCGCTGGACGCTGCATCGACTCCAGCGGGATCGGGGAACACACGGGCGAAGCGAAGTGCAGCAGGAGGCTTCGCCCAATCCTGGGGCGGGGGAAGATGCGCGCACCCGGACTGGCCCATAGGAACTCAAAGGAATTTCTGCGCAGAGCCTCCTGCTCTCAGCCTCGCTTCCTGCAGCTCCCTTTCCCTCTGGCGATAGTCACTAGTGTGCTTGAGTGACAGGCACCGGGAAGGAATAAACACCAGGACGCAAAAAGCACAGGGCTGGGCTGAGGCGGCTCCAGATTCAGACTCCTCTACCCGGGGTTGTGGGGTGGCGCGTGGGGTGGTGGGTGCTGGCTAAGCCCCAAGTCATGGACTCACAGTTTGTGTGACTGACTCTCCCCACTGGCCTAGAGGTCGTGGTTCACTGCTTGTGGCACCCCAGCTCTCCAAATGCTATCTGATGGAAGACTTGCTAGGTAAATGCTGGTTAAACAGGTAAACAAAGCAATAGCTGGAACCGGCTGGAGCATTCAACCTCCTCTGGGAGGTAGGCAGTCCTCTAACCATCACACAAGTACACACAGAGGTGCAGACCTCATCACAGGCCTGAGGAAGTCTGGAGGGTGTGGACAGAGGCATTCAAGGGGACCAGGGTTAGTCTGAGGGCTTCCTGGAGGAGGTGAGACTTAAGGATATGTATTAGGTGGAGGAGGTGGAGGAGGAGGAAGGCTTGGCTCAGGGAGGGGGAAGGACAGAGGAGACAAGGGCCCTGTGGTTGGAGAATGGAGGTGGGGAGATGTTGGGGTGATGAGGCTGGCCAGCTCCTGCAGGGGTGAGGAAGGAAGGGGTGGACTTCTGGGGAAGGGTGGGAGGTGCAGGGTGCCGAGTGGGTGCATCCTGTGCCATGAAAGTGATGTTTTTTTGGGGGGGATGGTGAAATCACTGATAAGTCATTTCCTGGGGTCTGTGAAGATCTTGCTTTGAGGAGTGCATCAGGTCAGCACCCTTCAGCCTGCTCCCTGCCCTGTCTCCCCAGGCCCTTACTCAGGCCCTGCCCTGGGTGTCCAGGTAAGCCAGGGTAAGCTGCCTTGGGCAGAGAAGGAAGTAGCTCCCGACCAGCTTAGCAATGCTTCCCCCGTTGGATCTTCTCATTTTCCCTGTCTGCATTTTTTTTTTTTTTTTTTTTGAGACAGAGTCTCGCTCTGTCACCCAGGCTGGAGTGCAATGGTGTGATCTTAGCTCACTGCAACCTCTGCCTCCTTGGTTCAAATGATTCTCCTGCCTCAGCCTCCGGAGTAGGTGGGTTTACAGGCATGGGCCACCATGCCCGGCTGATTTTGTATTTTTAGTATAGATGGGGTTTCTACATGTTGATCAGGGTGGTCTCAAACTCCCGACCTCAGGTGATCCTCCCACCTCGGCCTCCCAAAGTGCTGGGATTACAGGCGTGAGCCACCGTGCCTGTCGGATTTTAGTTTATTTTTTAAAGTGTTTCTGCTGAAATAACAATGCTCTGGGCTGTGAGGCGCAGGAAGAGCGGGTGGTGGAGGGAGAACTCCCCAACAATTCAGAGAGGTCATGTCCAGAGAGAGCGGTGACACGGGGTCCTCACCCTCCCCAAAGCTCAATCCAGTACTGTGAGTGTGGAGGGAGGCTTGCTGTGGGGCAGGGACTGGGTTGGAAGTGCTGCTCATGTTCCTAGGCTGGCTGAGAGCCAGATGGGCTGTGCAGCTTCAGCATGGATTCCAGGACTTCAGCACCGGGGACAGGGACGTTTCCAGAAGTTGCTGGCAGCACAGGAAACAGAAGCACCAGGAGGTCCCAGAGGCCCTGGGTAATTAGCTGGGAGCTCAGGGCTGGGCACACAGTGAGTTCTGAATGCCTTACAAGCTGCAGGAGACATTAGGAGTAAGGCCACCAAGGTGCTCTGCTGAAGAGAGTCTCCCCAGGTGTGGGCTGCCGGCTGTGAGCATCCCAGCCCTTCCTCCTTCCCCACTGACTGGGCTGCTCTCCAGAAGGGTCACCGTGACACAGGAAGTGTGAGTAGGGCCAGCACCCTCCTCCAGCCTACCCTCTGGCAGGAAGGGGTGGGAATGAGAGAATGTTCTGGCACCTGCACTTGCACTGGGGACAGCCCTATTTTGCTAGTTTGTTTTGTTTCGTTTTGTTTTTGATGGAGAGCGTATGTTAACTGGGTTCAGCCAGGAGATTTCAAAGGAGGAATTGTGATGCTGGGACTCCTATCAGGACAGGGGATGGTTCCAGCACCTTCTCTTGGGCAGGCCCACCTCAGCTTCCCCTCCCGCATACTCCTGCAGTTCTCCCTCCCCCAGCCCCATGCAGCTGCAGGAGGGTCTTGGTCTGAGACAGGTAAACACCTCCATTGTTGGCACATTCCGGGACAGAGAGAACCCAGGCACACAGGCTGCCGCCCTCAGACTTTATTCAAAGACCAGGAAGGGCCGGTGCAAGGAGGAGGAAGGGGCCCGTTGGGAGGCCCAGCGGGCGGAGGAACAGCTACCGAGCTCCGCGGCGGTGACGGTATTTGGAGGTTGGCACGGTGCTCACAGAAGCCGAAGACGTCAGGAGGCGTGCACCGCGGGGTGAACTCGGCGGGGAGGTGGGCGGCCCGGGGTCACCAGCGAGCGATGGCTTAGGAGCTGTGCAGAGGGTCGGTGCGGCCCAGGCCTGTGCCGCCTGCGCTACCTCCCGCAACCCGTCCTCTGCCCTGAAGGGAAGGCGCCATCTCGCCCCTCGACCCAGATCGCTCCCGGCCCGCCGCTCACCTTGAAGTTGACCGGGTCCACCCGAAGCTTGTGCGCGTGCAGGTCGCTCAGGGCGGACAGCGCGTTGGGCATGTCGTCCACGTGCGCCACGGCGTTTGGTCAGCGCGTCGGCCACCTTCTTGCCGTGGCCCTTAACCTGGGCAGAGCCGTGGCTCAGGTCGAAGTGCAGGAAGTAGGTCTTGGTGGTGGGGAAGGACAGGAACATCCTGCGGGGAGAAGCAGGAGTGAGGGGATGGGGTTTGGGTCCAGAACCCAGGACGGTTGAGGGTGGCCTGTGGGTCCGGGCGGGCCGAGGGAGCCCGGGTCGGAGCAGGGGAGGAGCCTCACCTCTCCAGGGCCTCCGCACCATACTCGCCAGCGTGCGCGCCCGACCTTACCCCAGGCGGCCTTGACGTTGGTCTTGTCGGCGGAGACAGCACCATGGTGGGTTCTCTCTGGTCTGTGGGGACCAGAAGAGTGCCGGGCCACGAGCGCGCCAGGGTTTATGCTTGGGGCGCGGGGCACGCCCGGCCGGGCGGCGCTCGTGGCTGGCGCGGAGCGGGGCGCGGCCTGGACGCAGGGGAGTCCCGGGCGGGGCCTGCGCGGGGCCGGGCGGGCTGGGCGGGGCCAGGACGTCTCCACCCTCCACCCGCCACTCCCACTCCCGCCCATCCCGCTCGCCCGGGAGCAGGGCGGTAGCGGCGCTGGCGGGGCCGGCCCGTTCGTGCACCCCGGGCGCAGCGCCGCCCTTTCCTTTCGGGCGCCGGAGCGTCCTAGCGGGGAGAAAAAGGAAGTCAGCCCGCACCCCGCCCCCCGGCCTGCCGCGCTGGACGCGCATCGACTCCAGCGGGATCGGGGAACACACGGGCGAGCGAGTGCGGCCGGGGCTTCCGCCCAATCCCTGGGGGCGGAGAGGGAATGGCGCACCCGGGACGCCTGGCCCATAGGAACTCAAAAGAGATTTCTGCGCAGAGCCCCTCCTGCTCTCCAGCCTCGCTTCCTACAGCTCCCTTTCCCTCTGGCGATAGTCACTGGTGTGCTTGAGTGACAGGCACCGGGGAAGGAACAAACACCAGGACGCAAAAAGCACGGGCTGGGCTGGGGGCGGCTCCAGATTTGAACTCCTCTACCCGGGGTTGGCCAGGGTGGCGCGTGGGGTGGTGGGTACTATGGCTAAGCCCCAAGTCCTTGGACTCACAGTTTGTGTGACTGGCTCTCCCACTGGCTAAGGGTCGTGGTTCACTGTGACACCCCCAGCTCTCCAAATGCTATCTGATGGAAGACATAGGTAAATACTGGTTGTAACAGATAAACAAACTTGGCTCTGGTGGGAAAGGACAGGGAAGGCCCCTGTGGTTTTATATGTAGGTGAGGGGATGTGGGGTGAGGGAAGGAGGGGTGGACTTGGCGAGGGGGAGGGTGGAGGTGCAGGAGTGCCAGTGGGTGCATCCCTACTCATGAAAAGTGATGGTTTGGCCAGGCTCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCAAGGCGGGTGGATCAGCTGAGGTCGGGAGTTCGAGACCAGGCTGACCAACATGGAGAAACCCCGTCTCTACTAAAAATACAAAAATTAGCTGGGTGTGGTGGTGCACGCCTGTAATCCCAGCTACTCGGGAGGTGGAGGCTGGAGAATCCTTTGAACCAGGGAGTTGGAGGTTGCAGTGAGCCGAGATCGTGCCACAGCACTCTAGCCTGACGACACAGCGAGACTCTGTCTCAAAAAAATAAATACATAAGGCCGGGCGCGGTGGCTCACGCCTGTAATCCCAGCACATTGGGAGGCTGAGGCGGGCAGATCACGAGGTCAGGAGATCGAGACCATCCTCGCTAATGCAGTGAAACTCTGTCTCTATTAAAAATAGAAAAAATTAGCCAGGCGTAGTGGAGGGCGCCTGTAGTCCCAGCTACTCAGCTACTCGGGAGGCTGAGGGAGGAGAATGGCGTGAACCTGCGAGGTGGAGCTTGCAGTGAGCCGAGATTGTGCCACTGCACTCCAGCCTGGGTGACAGAGCGAGACTCCGTCTCAAAAAAAAAAAAAAAATAATAATAATAATAAAAATAAAAAAAATAAAGTCATGATTTTTTGGGGGAGATCATGATGGAAACATAGTAATAATCAGTGAGACTGTGGAATGGCACAGAGCTGAATGAACTGGCTGAAGGATGCAGGAGACACTTCACTGAGAATAGGAAGTTGTACACGGGTATCTACAGTATGATGGTATTTTTGTCCAAAAGAGAGCCTGTGGCAGTAGTTGTAGATGTAGCTGTGTTCCCTCAGCATGGGATGGGGCCCATGGGCTGAGTTCAAACCCTTCATTTTGGTCTCTTCTGGGTGTGAGGAGACAGAAGAGAGACACTCTCCTACTTTAAAGTAACACAAAGGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGAAAGCATGGATAGATTTTAGTTTATTTATTTATTTATTTATTTATTTTTTTTGAGACAGAGTTTTGCTCTTGTTGCCCAGGCTGGAGTGCAATGGCGTGATCTCAGCTCACCTCAACCTCCACTCCCGGGTTCAAGCACTCCTGCCTCAGCCTTCCTAGTAGCTGGGATTACAGGCATGCGCCACCACACCCAGCTAATTTTGTATTTTTAGTATAGCTGGGGTTTCTCCATATTGGTCAGGCTGGTCTCAAACTCCCGACCTCAGGTGATCCTCTCGCCTCAGCCTCCCAAAGTACTGGGATTACAGGTGTGAGCCACCGTGTCGGATTTTAGTTTATTTTTTAAAGTGTTTCTGTTGAAATAACAATGCTCTGGGCTGTGAGGCACAGGAAAGGTGGGTGGTGGAGGGGAGAACTCCCCAACAATTCAGAGGTCCCATGTCCAGAAGAAAAGCGGTGACAGGGTCCTCTCCACCCTACCCAAAGCTCAATTCAGTACTGTGAGTGTGGAGGGGGGGGGGCTTGCTGGCGGGGCAGGGGGCTGGGCTGGAAGTGCTGCTCGTGTTCCTAGGCTGGCTGAAGGCAGATGGGCTGCCCAGCTTCCCAGCATGGATTCCAGGACTTCCAGCACCAGGGACAGGGACGTTTCAGAAGTTGCTGGCAGCACAGGAAGCAGAAGCAGGAGGCCGGTGGGGCAGTGGCTGGGGCTCAGGGCTGGGCACACAGTGAGTTCTGAATGCACCCTTCCTAGCTAGGAGACGTGCAGGGAGTAAGGCCACCATGAGTGCTGCACAGAGTCTCCCCAGGTGTGGGCTGCCGGCTGTGAGCATCCCAGCCCCTTCCTCCTTCCCCACTGACTGGGCTGCTCTCAGGGTCACCATGACACAGGAAGTGTGAGTAGGGCTAGCACCTCACTCAGCCACCTACCCTCTGGCAGGAAGGGGTGAGAGAATGAGAGAAATGTTCTGGCACCTGCACTTGCATTGGGGACCAGCCTGTTTGCTAGTTTGTTTTGTTTCGTTTTGTTTTGATGGAGAGCGTATGTTAACTGGGTTCAACCAAGGGAGATTTCAAAGGAGGGTGGCGATGCTGGGACTCCCTGAGTAGAACCCTGGATGCCTCTCTCCCTCCCCAGCTTCTGAGCCACTGCCTGCAGGCCTGGCCTCTCAGGACAGAGGATGGTTCAGCACCTTCTCTTGAGCAGGCCCACCTCCAGCTTCCCCTCCATACTCCCTGCAGTTCTCCCTCCCCAGCCCATGCAGCTGCAGGAGGTCCTTGGTCTGAGACAGGTAAACACCTCCACTGTTGGCACATTCCAGGACAGAGGGAACCCAGGCACACACAGGCTGCTGCCTACTCGGACTTCATTCCAAAGATCAGGAAGTGCTGGGGCAGGGAGGGAGGAGGGCCTGGTGGGAGGCCCAGCAGGCAGGGGGCACTGGTGGCCACCAATGCTCTAGCTTAATGGTGTTTTAAGGTCAACAGCACTCACAAAAGCAAATAATTTGTCCATGGAGGCACCACATGGGTAAATGGGAAGGGTGTGGCAGGCCTAGGTCACCAGCAGGCATTTACCCAGGAACTGTGGAGAGAGAGTCAGGGGATGCAGATACTACTGCCACCTGTGCTCCACCTCCCTGAACGCGTGGGAGGGTCTGTCCTATGAGGAAGGCCAGCTGCCCCCACACCCTAGATTCCTCCCAGCCCAAAGCTCTCCCTGATGTTGCCTGGGTCCACCCACAGCTCGTGGACATGCAGCTTCCTCAGCTCAGACACATCATTGGGCATGTCATCTAAGTGGACCACGGCATTGGTCAGCGTGTCGGCCACCTACTTGCATTGTTTCTTGATCTGGGCTGAGCCATAACTCAGGGCCAAGTGAGAGGAAGGTCAGGAACATGCTTCAGGGGGAAGCAAAGAGTGAGGGGCTGGGGGCTGGGTCCAGGGTAAGAATAGTTAGGGGTGGCCAGTGGAGCCCTGGAGTGGGCGTTGGACCCAAGGTGTGAGCAGGGAGAGGAGGGTTCTTGCCTCTCCAGGGCCTCCGTGGCATAGCCAGCAGTGTGGTTGCCAACTTTCTCCCAGACCGCCTTGGTGTTAGCCTTGTCCTCAGGAGACAGCACCACGGCAGCTTATTTCTGAGTCTGTGTAAGCTCAAAAGTGAAAGGTGGGGGAGAGTAGGAGCTCTTATGCCTGGGATGCAGGGTCATTGGCTGGCAAGTTGCCAGGAGGTGCCTGGCACAGGAAAGGGGAGGAGGAGTGAGTGAGAGTAGACAGGTTGCCCAAGGCTGCCCTGCGTCTCTTAGTCCTGCCCAAACAGAGACCAGCTCTCCACATGGTAGCCTTCCCGCTGGCATTGTCCTGAGCTTCCACACAGGGTCCTGAGCAGCTTATACCCATAAGGTATCCACATATGGCCCAGACCCAGGCAGCTTCCTGGATGGGGAGAAGAGCACAGTTTTTCCAAAGTGAAACATAAGCTGCTTTCACACCAAACTCTGCCTGGGACAAAAAGCTGTTAATTATCCATGTAAACAGTTGCGGGGGCTGGGTGTGGTGGCTCACGCCTGTCATCCCAGCACTTTGGGAGACCGAGGAGGGCAGATCACTTGAGGTCAGGAGTTTGAGACCAGCCTGATCAACATGGTGAAACTGTGTCTCTATTAAAAATACAAAAATTAGCCAGGAATGGCGGCCCATGCCTGTAATCCCAGCTACACGGGAGACAGAGGCAGGAGAATCGCTTGAACCCAGGAGGTTGAGGCTGCAGTGAGCCAAAACTGTGCCACTGCACTCCAGCCGGGGTGTCAGAGCAAGGGCCTATCTCAAAAACAAAAACAAAGCAAAAAAAAAAAAAACCCCCAAAAACAAAAAAACAGTTGGGGGATGTAGATAACGTGGGGCCAAGTTCTCCAAACTACCGGGCTGGGGTCCCAGGCAGTGCAGGGCAACAGGTCACTGGCACTGACTGCTGGGGGGAGTTTGGCCGACTCACAGTCAGGGCTCCACTTCCAGCAGCTCTCTTCAGACCTCGTTGCTGACTGCCAGGTCAGACCTGGGTGCCAGGCCTGGTCCAGTGCTCCCCACACAGACCCAGGATGGCCAACCCAGAGCCAACCCTATTCAGGGTTTGGAGACCTCCAAGGTCCCACTCTTTCCTCAGCCCCTATTCTTTGTCTGAAAAGCCTGGGGTGGAGCCCTAAGATTCCCTTCCCAGTCCTGGGGAGGGAGAAGAGGGAACCTGGGGAAGGGAGTGCCTTGGCCTGGGCACCTGCAGAGATTGAGACACGGCAGCTGGGGCTGGGTGGAGGGGTGGGGTGGCAACGACCCTGCAGGCACCCACCCTGTGTTATGATTTCAAAAATTATTTCAAAAATCTCGTGAGTTGGACGTTGTTTTTCTATTCGATAGGTAAAGGCCTGAGAGCTCTGGGACACGCCCAGGGTCGCACAGCTTGTCATTGAGCAGGCTGGAAGACTAAACCTGCTCTGTCTGCCAAATTGTGGCCACTGAGGACTTCCTCCCGACCTCATCCC";
	let z = b"GGGATGAGGTCGGGAGGAAGTCCTCAGTGGCCACAATTTGGCAGACAGAGCAGGTTTAGTCTTCCAGCCTGCTCAATGACAAGCTGTGCGACCCTGGGCGTGTCCCAGAGCTCTCAGGCCTTTACCTATCGAATAGAAAAACAACGTCCAACTCACGAGATTTTTGAAATAATTTTTGAAATCATAACACAGGGTGGGTGCCTGCAGGGTCGTTGCCACCCCACCCCTCCACCCAGCCCCAGCTGCCGTGTCTCAATCTCTGCAGGTGCCCAGGCCAAGGCACTCCCTTCCCCAGGTTCCCTCTTCTCCCTCCCCAGGACTGGGAAGGGAATCTTAGGGCTCCACCCCAGGCTTTTCAGACAAAGAATAGGGGCTGAGGAAAGAGTGGGACCTTGGAGGTCTCCAAACCCTGAATAGGGTTGGCTCTGGGTTGGCCATCCTGGGTCTGTGTGGGGAGCACTGGACCAGGCCTGGCACCCAGGTCTGACCTGGCAGTCAGCAACGAGGTCTGAAGAGAGCTGCTGGAAGTGGAGCCCTGACTGTGAGTCGGCCAAACTCCCCCCAGCAGTCAGTGCCAGTGACCTGTTGCCCTGCACTGCCTGGGACCCCAGCCCGGTAGTTTGGAGAACTTGGCCCCACGTTATCTACATCCCCCAACTGTTTTTTTGTTTTTGGGGGTTTTTTTTTTTTTTGCTTTGTTTTTGTTTTTGAGATAGGCCCTTGCTCTGACACCCCGGCTGGAGTGCAGTGGCACAGTTTTGGCTCACTGCAGCCTCAACCTCCTGGGTTCAAGCGATTCTCCTGCCTCTGTCTCCCGTGTAGCTGGGATTACAGGCATGGGCCGCCATTCCTGGCTAATTTTTGTATTTTTAATAGAGACACAGTTTCACCATGTTGATCAGGCTGGTCTCAAACTCCTGACCTCAAGTGATCTGCCCTCCTCGGTCTCCCAAAGTGCTGGGATGACAGGCGTGAGCCACCACACCCAGCCCCCGCAACTGTTTACATGGATAATTAACAGCTTTTTGTCCCAGGCAGAGTTTGGTGTGAAAGCAGCTTATGTTTCACTTTGGAAAAACTGTGCTCTTCTCCCCATCCAGGAAGCTGCCTGGGTCTGGGCCATATGTGGATACCTTATGGGTATAAGCTGCTCAGGACCCTGTGTGGAAGCTCAGGACAATGCCAGCGGGAAGGCTACCATGTGGAGAGCTGGTCTCTGTTTGGGCAGGACTAAGAGACGCAGGGCAGCCTTGGGCAACCTGTCTACTCTCACTCACTCCTCCTCCCCTTTCCTGTGCCAGGCACCTCCTGGCAACTTGCCAGCCAATGACCCTGCATCCCAGGCATAAGAGCTCCTACTCTCCCCCACCTTTCACTTTTGAGCTTACACAGACTCAGAAATAAGCTGCCGTGGTGCTGTCTCCTGAGGACAAGGCTAACACCAAGGCGGTCTGGGAGAAAGTTGGCAACCACACTGCTGGCTATGCCACGGAGGCCCTGGAGAGGCAAGAACCCTCCTCTCCCTGCTCACACCTTGGGTCCAACGCCCACTCCAGGGCTCCACTGGCCACCCCTAACTATTCTTACCCTGGACCCAGCCCCCAGCCCCTCACTCTTTGCTTCCCCCTGAAGCATGTTCCTGACCTTCCTCTCACTTGGCCCTGAGTTATGGCTCAGCCCAGATCAAGAAACAATGCAAGTAGGTGGCCGACACGCTGACCAATGCCGTGGTCCACTTAGATGACATGCCCAATGATGTGTCTGAGCTGAGGAAGCTGCATGTCCACGAGCTGTGGGTGGACCCAGGCAACATCAGGGAGAGCTTTGGGCTGGGAGGAATCTAGGGTGTGGGGGCAGCTGGCCTTCCTCATAGGACAGACCCTCCCACGCGTTCAGGGAGGTGGAGCACAGGTGGCAGTAGTATCTGCATCCCCTGACTCTCTCTCCACAGTTCCTGGGTAAATGCCTGCTGGTGACCTAGGCCTGCCACACCCTTCCCAGTTTACCCATGTGGTGCCTCCATGGACAAATTATTTGCTTTTGTAGGTGCTGTGTTGACTAAAAACACCATTAAGCTAGAGCATTGGTGGTCATGCCCTGCCTGCTGGGCCTCCTGGCCCTCTTTCTCCCCCTCCCTGCCCAGCACTTCCTGATCTTTAGATGAAGTCCGAGTAGGCAGCAGCCTGTGTGTGCCTGGGTTCTCTGTCCCGGAATGTGCCAACAGTGGAGGTGTTTACCTGTCTCAGATTGGGACCTCTCTGCAGCTGCATGGGGGCTGGGGAGGGGGAGAACTGGGGAGTATGGGAAGCTGAGGTGGGCCTGCTCAAGGAAAGGTGCTGAACCATCCCCACCATCCTGAGAGGTGCCAGGCCTGCAGGCAGTGGCTCAGAAGCTGGGAGGAGAGAGGCATCCAGGGTTCTACTCAGGAGTCGCCACCCTTCTTTGAAATCTCCCTGGTTGAACCCAGTTAACATACGCTCTCCATCAAAACAAAACGAAACAAAACAAACTAGCAAAATAGGCTGTCCCCAATGCAAGTGCAGGTGCCAGACATTTCTCTCATTCATCACCCCTTCACTGCCAGAGGGTAGGTGGCTGGAGTGAGGGTGCTGAGCGCCCCTACCTCACACTTCCTGTGTCATGGTACCCTCCTGAGAGCAGCCCAGTCAGTGGGAAGGAGGAAGGGGCTGGGATGCTCACAGCCGGCAGCCACACCTGGGGAGACTCTTCAGCAGAGCACCTTGCGGCCTTACTCCTCGTCTCCTGCAGTTTGTAAGGTGCATTCAGAACTCACTGTGTGCCCAGCCTGAGCTCCCAGCTAATTGCCCCACCCAGGGCCTCTGGGACCCTCCTGGTGCTTCTGCTTCCTGTGCTGCCAGCAACTTCTGGAAACGTCCCTGTCCCCGGTGCTGAAGTCTGGAATCCATGCTGGGAAGTTGCACAGCCCATCTGGCTCTCAGCCAGCCTAGGAACACGAGCAGCACTTCCAGCCCAGCCCCTGCCCCAGCAAGCCTCCCCCTCCACACTCACAGTACTGAATTGAGCTTTGGGTAGGGTGGAGAGGACCCTGTCACACACTTTTCTTCTGGACATGGACCTCTCTGAATTGTTGGGGGAGTTCCCTCCCCCTCTCCACCACCCACTCTTCCTGTGCCTCACAGCCAGAGCATTGTTATTTCAACAGAAACACTTTAAAAAATAAACTAAAATCCGACAGGCACGGTGGCTCACACCTGTAATCCCAGTACTTTGGGAGGCTGAGGCAGAAGGATCACCTGAGGTCGGGAGTTTGAGACCAGCCTGACCAATATGGAGAAACCCCAGTTATACTAAAAATACAAAATTAGCTGGGTGTGGTGGCGCATGCCTGTAATCCTAGCTACTAGGGAGGCTGAGGCAGGAGAATCGCTTGAACCCGGGAGGTGGAGGTTGAGGTGAGCTGAGATCACGCCATTGCACTCCAGCCTGGGCAACAAGAGCAAAACTCCGTCTCAAAAAATAAATAAATAAATAAATAAATAAACTAAAATCTATCCATGCTTTCACACACACACACACACACACACACACACACCCTTTTTGTGTTACTTAAAGTAGGAGAGTGTCTCCTCTTTCCTGTCTCCTCACACCCACCCCCAGAAGAGACCAAAATGAAGGGTTTGGAACTCAGCCCATGGGCCCCACCATCCCATGCTGAGGGGAACACAACTACATCTACAACTCGCACAGGCTCTCTTTTTGGACAAAAATACCATCATACTGTAGATACCCGTGTACAACTTCTATTCTCAGTGAAGTGTCTCCCTACATCCCTTTCAGCCAGTTCATTCAGCTCTGTGCCATTCCACAGTCTCCTGATTATTTCTTTTTTTTTTTTTTTTTTTTGAGACGGAGTCTCGCTCTGTCACCCAGGCTGGAGTGCAGTGGCACAATCTCGGCTCACTGCAAGCTCCACCTCGCAGGTTCACGCCATTCTCCTCCCTCAGCCTCCCGAGTAGCTGAGTAGCTGGGACTACAGGCGCCCCCCACTACGCCTGGCTAATTTTTTCTATTTTTAGTAGAGACAGAGTTTCACTGCATTAGCGAGGATGGTCCGATCTCCTGACCTCGCATCTGCCCGCCTCAGCCTCCCAATGTGCTGGGGTGGGATTACAGGCGTGAGCCACCGCGCCCGGCCTTGTATTTATTTTTTTGAGACAGAGTCTCGCTGTGTCGTCAGGCTAGAGTGCTGTGGCACGATCTCGGCTCACTGCAACCTCCAACTCCCTGGTTCAAAGGATTCTCCAGCCTCCACCTCCCGAGTAGCTGGGATTACAGGCGTGCACCACCACACCCAGCTAATTTTTGTATTTTTAGTAGAGACGGGGTTTCTCCATGTTGGTCAGCCTGGTCTCGAACTCCCGACCTCAGCTGATCCACCCGCCTTGGCCTCCCAAAGTGCTGGGATTACAGGCGTGAGCCACCGAGCCTGGCCAAACCATCACTGCTTCATGGCTGGGGATGCACCCACTGGCACTCCTGCACCTCCCCACCCTCCCCTCGCCAAGTCCACCCCTTCCTTCCTCACCCCACATCCCCTCACCTACATTCTGCAACCACAGGGGCCTTCTCTCCCTGTCCTTTCCCTACCAGAGCCAAGTTTGTTTATCTGTTTACAACCAGTATTTACCTAGCAAGTCTTCCATCAGATAGCATTTGGAGAGCTGGGGGTGTCTGGTAGACCACGACCTCTAAGGCCAGTGGGAGAGTCCAGTCCCCACACAAACTGTGAGTCCATGACTTGGGGCTTAGCCAGTACCCACCACCCCGCGCCACCCCACAACCCGGGTAGAAGAATTTGAATCTGGGCCGCCCCCAGCCCAGCCCCGTGCTTTTTGCGTCCTGGTGTTCCTTCCCGAGTGCCTGTCACTCCAAGCACGCAAAAGTGACTATCGCCAGAGGAAAGGGAGCTGCAAGGAAGCGAGGCTGGAGAGCAGGAGGGGGCTCTGCGCGAAATTCTTTTGAGTTCACTATGGGCCAGGGCGTCCGGGTGCGCGCATTCCTCTCCGCCCCAGGATTGGGCGAAGCCTCCCGGCTCGCACTCGCTCGCCCGTGTTCCCCGATCCCGCTGGAGTCGATGCGCGTCCAGCGCGTGCCAGGCGGGGCGGGGGTGCGGGCTTGTTGACTTTCTCCCTCGCTAGGGACGCTCCGGCGCCCGAAAGGAAGGGTGGCGCTGCGCTCGGGGTGCACGGGCCGACAGCGCCCGGCCCAACGGGCCCGGCCCCGCCAGCACCTGCCGCCCTGCTCCCGGGCGAGCGGGATGGGCGGGAGTGGAGTGGCGGGTGGAGGGAGTGGAGACGTCCTGGCCCCCGCCCCGCGTGCACCCCAGGGAGGCCGGGCCGCCGCCCGGCCCCGCGCGGGCCCCGCCCGGGACTCCCCTGCGGTCCGGGGCACCCGGGCTCCGCGCCAGCCCAATGAGCGCCGCCCGGCCGGGCGTGCCCCCGCGCCCCAAGCATAAACCCTGGCGCGCTCGCGGCCCGGCGCTCTCTGGTCCCCACAGACTCCAGAGAGAACCCACCATGGTGCTGTCTCTGCCGACAAGACCAACGTCAAGGCCGCCTGGGGTAAGGTCGGCGCGCACGCTGGCGAGTATGGTGCGGAGGCCCTGGAGAGGTGAGGCTCCCTGCCTGCTCCGACCCGGGCTCCCCGCCCGGACCCACAGGCCACCCTCAGCCGTCCTGGCCCCGGACCCAAACCCCACCCTCGCTCTGCTTCTCCCACAGGATGTTCCTGTCCTTCCCCACCACCAAGACTACTTCCGCACTTCGACCTGAGCCACGGCTCTGCCCAGGTTAAGGGCCACGGCAAGAAGGTGGCCGGCCTTGTGACCAGCACCGTGGCTACGTGGACGACATGCCCAACGCGCTGTCCGCCCTGAGCGACCTGCACGCGCACAAGCTTCGGGTGGACCCGAGTCAGCTTCAAGGTGAGCGGCGGAGCCCGAGGCGATCTGGGTCGAGGGGCGAGGGCGCCTTCCTCTCAGGGCAAGAGGATCACATGAGTTGCGGGAGGTGTAGCGCAGGCGGCGGCTGCGGGCCTGGGCCCGCACTGACCCTCTTCTCTGCACAGCTCCTAAGCCTGCCTGCTGGTGACCCTGGCCGCCCACCTCCCCGCCGAGTTCACCCCTGCGGTGCACGCCTCCTAAGTTCCTGGCTTCTGTGAGCACCGTGCTGACCTCCCAAATACCGTCAAGCTGGAGCCTCAGTAGCCGTTCCTCCTGCCCGCTGGGCCTCCAAGCAGGCCCTCCTCCCCTCCTTGCACCGGCCCTTCCTGGTCTTTGAATAAAGTCTGAGTGGGCGGCAGCCTGTGTGTGCCTGGGTTCTCTCTGTCCCGGAATGTGCCAACAATGGAGGTGTTTACCCCTGTCTCAGACCAAGGACCTCTCTGCAGCTGCATGGTGGGGAGGGAGGGCACTGCAGGGGGTAGCAGAGCTGAGTGGGCCTGCTCCAGAGAAGGTGCTGAACCATCCCCTGTCCTGAGGTGCCAGGCCTGCAGGCAGTGGCTCCAGAAGCTGGGGAGGAGAGAGGCATCCAGGGTTCTACTCAGGGAGTCCCAGCATCGCCACCCTCCTTTGAAATCTCCCTGGTTGAACCCAGTTAACATCGCTCTCCATCAAAACAAAACGAAACAAAACAAACTAGCAAAATAGGCTGTCCCCAGTGCAAGTGCAGGTGCCAGAACATTTCTCTCTATTTCCACCTTCCTGCCAGAGGGTAGGTGGCTGGAGTGAGGGTGCTGGCCCTGCTCACACTTCCTGTGTCACAGTGACCCTCTGAGAGCAGCCCAGTCAGTGGGAAGGAGGAGGGGCTGGGATGCTCACAGCCGGCAGCCCACACCTGGAGACTCTTCAGCAGAGCACCTTGCAGCCTTACTCCACCTGCGTCTCCCTGCAGTTTGTAAGGTGCATTCAGAACTCCTGTGTGCACCCCAGCCCTGAGCTCCCAGCTAATTGCCCCCAGGGCCTCTGGGACCTCTGAGTGCTTCTGGCCTAACTTGCCAGCAACTTCTGGAAACGTCCCTGTCCCCGGTGCTGAAGTCACCAAGTCCATGCTGGGAAGTTGCAGCCCATCTGGCTCTCAGCCAGCCTAGGAACATGAGCAGCACTTCCAACCCAGTCCCTGCCCCACAGCAAGCCTCCCCCTCCACACTCACAGTACTGGATTGAGCTTTGGGAGGGTGGGAGAGGACCCTGTCACCGCTTTCCTTCTGGACATGGACCTCTCTGAATTGATGCCTCCCCCTCTCACCCACCCGCTCTTCCTGCGCCTCACAGCCCAGAGCATTGTTATTTCAGCAGAAACACTTTAAAAAATAAACTAAAATCCGGCTGGGCACGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGGTGGGAGGATCACCTGAGGTCGGGAGTTTGAGACCACCCTGATCAACATGTAGAAACCCCATCTATACTAAAAATACAAAATCAGCCGGGCATGGTGGCCCATGCCTGTAAACCCCACCTACTCCGGAGGCTGAGGCAGGAGAATCATTTTAACCAAGGAGGCAGAGGTTGCAGTGAGCTAAGATCACACCATTGCACTCCAGCCTGGAAAACAACAGCGAAACTCCGCCTCAAAAAAAAAAAAAGCCCCCACATCTTATCTTTTTTTTTCCTTTCAGACTGTGGGCAGAGTCAGAAAGAGGGTGGCAGACAGGGAGGAAATGAGAAGATCCAACGGGGGAAGCATTGCTAAGCTGGTCCAGGCTACTTCCTTCTCTGCCCAAGGCAGCTTACCCCCAGCTGCTCCTGGACACCAGGGGCAGGGCCTGAGTAAGGGCCTGGGAGACAGGAGCAGAGAGCAGGCTGAAGGGTGCTGACTGATGCACTCCTCAAAGCAAGATCTTCTGCCAGACCCCCAGGAAATGACTTATCAGTGATTTCTCAGGCTGTTTTCTCCTCAGTACCATCCCCCCAAAACATCACTTTTCATGCACAGGGATGCACCCACTGCTCCTGCACCTCCCACCCTTCCCCAGAAGTCCCACCCCTTCCTTCCTCACCCTGCAGGAGCTGGCCAGCCTCACCCCAACATCTCCCCACCTCCATTCTCCAACCACAGGGCCCTTGTCTCACTCTGTCCTTTCCCCTCCCCGAGAACCAAGCCTCCTCCCTCCTCCACACCCTCCTCCACCTAATACATATCCTTAAGTCTCCTCCTCCAGGAAGCCCTCAGACTAACCCTGGTCCCCTTGAATGCCTCGTCCACACCTCCCAAGACTTCCTCCAGGGCCTGCTTGGTAAATTCTGCACCTCTGTGTGTACTTGTGTGATGGTTAGGGGACTGCCTACACTCCAGAGGAGGTTGAATGCTCAACCGGTTCCAGCTATTGCTTTGTTTACCTGTTTAACCATTACCCACACTTCCTGCAACCCCGGGTAGAGGAGTCTGAATCTGGAGCCGCCCCCAGCCCAGCCCGTGCTTTTTGCGTCCTGGTGTTTATTCCTTCCCGGTGCCTGTCACTCAAGCACACTAGTGACTATCGCCAGAGGAAAGGGAGCTGCAGGAAGCGTGGTTCCTATGGGCGGGGCGTCCGGGTGCGCATTCCTCTCCGCCCCAGGATTGGGCGAAGCCTCCGGCTCGCCGCTCAGCTCGCCCGTGTGTTCCTACGATCCCGCTGGAGTCGATGCGCGTCCAGCGCGTGCCGGGCCGGGGCGGGGTGCGGGCTGACTTTCTCCCTCGCTAGGGACGCTCCGGCGCCCGAAAGAAAGGAAAGGGTGGCGCTGCGCTCCGGGTGCGAGCCGACAGCGCCCGACCCCAACGGGCCGGCCCCGCCAGCGCCGCTACCGCCCTGCTCCCGGGCGAGCGGGATGGCGGGAGTGGAAGTGGCAGGTGGAGGGTGGGAGACGTCCTGGCCCCCGCCCCGCGTGCACCCCCAGAGGCCGAGCCCGCCGCCCGGCCCCGCGCGGGCCCCGCCGGGACTCGGTCAGAGCCCCGCGCCCCGGGCTCCGCGCCGGCCAATGAGCGCCGCCCGGCCGGGCGTGCCCCCGCGCCCAAACGCAAACCCTGGCGCGCTCGCGGCCGGCACTCTTCTGGTCCCCACAGACTCAGAGAGAACCCACCATGGTGCTGTCTCCTGCCGACAAGACCAACGTCAAGGCCGCCACAGAGCAAGGTCGGCGCGCACGCTGGCGAGTATGGTGCAGGCCCTGGAGAGGTGAGGCTCCCTCACCCTGCTCCGACCCGGGCTCCTCGCCCGCCCGGACCCCAGGCCACCCTCAACCGTCCTGGCCCCGGACCCAAACCCACCCCACCCCTCACTCCACCTGCTTCTCCCCGCAGGATGTTCCTGTCCTTCCCCACCACCAAGACCTACTTCCCGCACTTCGACCTGAGCCACGGCTCTGCCAGGTTAAGGGCCACGGCAAGAAGGTGGCCGACGCGCTGACCAACGCCGTGGCTACGTGGACAGCATGCCCAACACGCTGTCCGCCTGAGCGACCTGCACGCGCACAAGCTTCGGGTGGACCCGGTCAACTTCAAGGTGAGCGGCGGGCCGGGAGCGATCTGGGTCGAGGGGCGAGATGGCGCCTTCCTCGCAGGGCAGAGGATCACGCGCGGGTTGCGGGAGGTGTAGCGCAGGCGGCGGCTGCGGGCCTCGGCCCCACTGACCCTCTTCTCTGCACAGCTCCTAAGCCACTGCCTGCTGGTGACCCTGGCCGCCCACCTCCCCGCCACGAGTTCACCCTGCGGTGCACCTCCCTGGACAAGTTCCTGGCTTCTGTGAGCACCGTGCTGACCTCCAAATACCGTTAAGCTGGAGCCTCGGTGGCCATGCTTCTTGCCCCTTGGGCCTCCCCCCAGCCCCTCCTCCTTCCTGCACCCGTACCCCGTGGTCTTTGAATAAAGTCTGAGTGGGCGGCAGCCTGTGTGTGCCTGGGTTTTCCCTCAGCAAGCGTGCCAGGCATGGGCGTGGACAGCAGCTGGGACACATGGCTAGAGACCTCTCTGCAGCTGGATAGGGTAGAGGAAAGGCAGGGGCGGGAGGAGGGATGGAGGAGGGAAAGTGGAGCCACCGCGAAGTCCAGCTGGAAAAACGCTGGACCTAAGAGTGCTTTGAGGATGCATTTGCTCTTTCCCGAGTTTATTCCCAGACTTTTCAGATACAGGTTTGCTGAAATAATGAATTTATCCATCTTTACGTTTCTGGGCACTCTTGTGCCAGGGGAGCAGCTGGCTTTCTGCCTGGGACGTCACTGGTTTCCCAGGGTCCTCCCACATATGGGTGGTGGGTAGGTCAGAGAGTCCTCAGCATGGCTGCATTGATCCCCCATCGTTCCCACTAGTCTCCGTAAAACCTCCCAGATACAGGCACAGTCTAGATGAAATCAGAGGGTGCGGGGTGCAGCACAGGCCCCAGGCAATTCAATAGGGGCTCTACTTTCACCCCCAGGTCACCCAGAATGTCACAGACATGACGCCCTGGGGCTGTCAAGATCAGGCGTTTGTCTCAGGGCCCAGCTCAGGGCCCAGCTCAGCACCCACTCAGCTCCCTGAGGCTGGGGAGCCTGTCCCATTGCGACTGGAGAGGAGAGCGGGGCCACAGGGAGCCTGGCTAGAAGGTCCCTTCTCCCTGGTGTGTGTTTTCTCTCTGCTGAGCAGGCTTGCAGTGCCTGGGGTATCAGAGGTGAGGTTCAGAGCTGGTAACAAAGCCTGGCCCTCAACTGATAGGAATATCTTGCATTCCCTGAGCCCATGAATCTACCCTTGGTAAGCACACCATATGGCAGGCCCTCTGCTGCGTTTGTGATGTCCTTCCCGCAGCCTGTGGGTACAGTATCAACTGTCAGGAAGACAGTGTCTTCGTTATTTCATCCGGAAGAATGGAGGTCTGACCTAAAGGTAGAAATATGTCAAATGTACAGCAGAGGCTGGTTGGGTGTAGCGCTTTTTACAATTAATTGATCAGAACCAGTTATAAATTTATCATTTCCTTCTCCACTCCTGCTGCTTCAGTTGACTAAGCTAAGAAAAAATTATAAAAATTGGCCGGGCGCGGTGGCTCACACCTGTAATTGCAGCACTTTGCCAGGCTTAGGCAGGTGGATCACCTGAAGTCAGGGGTTCGAGACCAGCCTGGCCAACATAGTGAAACCCTGTCTCTACTAAAGACAAAAATTGTCCAGGTGTGCCACCTCTGTAAACCTGGCACTTTGGGAGGCGGAGGTTGTAGTGAGTCAAGATCGCGCCATCGCACTCCAGCTTGGGCAACAAGAGCGAAACTCTGTCTCAAAAAAAAATTTAATCTAATTTAATTTAATTTAAAAATTAGCACGGTGGTTGGGCACAGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAAGCCAAGGTGGGCAGATCACAAGGTCAGGAATTCGAGACCAGCCTGGCCAATATGGGGAAACCCCATCTCTACTAAAAATACAAAAAATTAGCCGGGTGTGGTGGCGCACGCCTGTAATCCCAGCTACTCGGGAGGTTGAGGTAGGAGAATCACTTGAACCCAGGAGGCAGAGGTTGCAGTGACCCGAGATCACACCATTGCACTCTAGCCTGGGCAACAAGAGCAAAACTCCATCTCAAAAAAAATTATAAAAATTATACATCAGTAGATGAATGGGTAAACAAAATGTGGTGGTCTATACACACAATGGAATATTATTTGGCCACAAAAAGAAATGAAGCACTGATAGGATGTAGCTGCACCCTGAAAATATTTGACAAGTAAAAGAAGCCGGACACCAAAGGTCACAAACTGCATGACCCCATCTATATGCAATATCCGCTACAGCCAAATCCATAGGGACCAAAAGCGGATTAGTGGCTGCCGGGGCCAGAGTTACTGTTAATGAGTACCGAGGTGGCGTTTGGGATGATGAAAAAGTTCTGACCTAGATAGTGGTGATGGCTGCATAACACTAAGTGTTCTTAATATCACCAAATTTTATACCTGAAAAATGGCTACAATGGTAATTTATGTCTATTTTATCACCTTTTTTAAAACAAAAAAGATATAAGGGGTACAGCAGAGTGAGTGCTGCATATGCATTTACTATTATTCTTGGGTTACATCCCAGGTACTCAATAAATGTTCACTGCCCTGAAGAAACACCTGCTACGAGTCAGGCACCTCACAGTTGTTATCCGTTTAATTCTCACAATCTGAGAAGAAACTGTCACCCTCATTTTATATAATAAATGAGAAAACAGACTCGGGCAAGTGTCACAATAGAATCAAGAGGCAGAATAAACTGACTTCCAATGCCAAATCCATGCCGAAATTCAGTGCTATAATAATGTACATGGCCGGGCGCGGTGGTTCACGCCTGTAATCCCAGAACTTTGGGAGGCTGAGGCGGGAGGATCACCTGAGGTCGGGAGTTTGAGATCAGCCTAACACGGTGAAACCCTGTCTCTACTAAAAATACAAAATTGGCATGGTGGCATGCACCTGTGATCCCAGTTACTCGGGAGGCTGAGGCAGGAGAATCGTTTGAACCCGGGAGGCGGAGGTTGCAGTGAGCCGGAATGGCGCCACTGCACTCACCGCACCCGGCCAATTTTTGTGTTTTTAGTAGAGACTAAATACCATATAGTGAACACCTAAGACGGGGGGCCTTGGATCCAGGGCGATTCAGAGGGCCCCGGTCGGAGCTGTCGGAGATTGAGCGCGCGCGGTCCCGGGATCTCCGACGAGGCCCTGGACCCCCGGGCGGCGAAGCTGCGGCGCGGCGCCCCCTGGAGGCCGCGGGACCCCTGGCCGGTCCGCGCAGGCGCAGCGGGGTCGCAGGGCGCGGCGGGTTCCAGCGCGGGGATGGCGCTGTCCGCGGAGGACCGGGCGCTGGTGCGCGCCCTGTGGAAGAAGCTGGGCAGCAACGTCGGCGTCTACACGACAGAGGCCCTGGAAAGGTGCGGCAGGCTGGGCGCCCCCGCCCCCAGGGGCCCTCCCTCCCCAAGCCCCCCGGACGCGCCTCACCCACGTTCCTCTCGCAGGACCTTCCTGGCTTTCCCCGCCACGAAGACCTACTTCTCCCACCTGGACCTGAGCCCCGGCTCCTCACAAGTCAGAGCCCACGGCCAGAAGGTGGCGGACGCGCTGAGCCTCGCCGTGGAGCGCCTGGACGACCTACCCCACGCGCTGTCCGCGCTGAGCCACCT";
	let scoring = Scoring::new(-1, 0, |a: u8, b: u8| if a == b { 3 } else { -4 });
	let mut aligner = Aligner::new(scoring, x);
    let sequence: Vec<&[u8]> = vec![y,z,p];
	// let g = aligner.graph().map(|_, n| (*n) as char, |_, e| *e);
	// println!("Dot graph:\n{:?}", Dot::new(&g));
	// println!("start poa");
	// let coy = aligner.global(y).add_to_graph_and_correct();
	// let coy_seq = get_node_weights(&aligner.poa.graph, &coy);
	// println!("{:?}", aligner.alignment());
	// println!("node_indices: {:?}", coy);
	// println!("node_indices_id: {:?}", coy_seq);
	// println!("node_indices_seq: {:?}", String::from_utf8_lossy(&coy_seq));
	// let coy = aligner.global(p).add_to_graph();
	// println!("{:?}", String::from_utf8_lossy(&coy));
	// println!("{:?}", aligner.alignment());
	for seq in sequence.clone(){
		aligner.global(seq).add_to_graph();
		println!("{:?}", aligner.alignment());
	}
	let consensus = aligner.consensus();
	println!("{:?}", consensus);
	println!("{:?}", String::from_utf8_lossy(&consensus));
	let g: Graph<char, i32, Directed, usize> = aligner.graph().map(|_, n| (*n) as char, |_, e| *e);
    println!("Dot graph:\n{:?}", Dot::new(&g));
	println!("{:?}", aligner.alignment());
}
