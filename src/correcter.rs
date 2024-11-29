use petgraph::{Directed, Graph};
pub type POAGraph = Graph<u8, i32, Directed, usize>;
use petgraph::graph::NodeIndex;
use petgraph::visit::EdgeRef;
pub type Indices = Vec<NodeIndex<usize>>;
use log::{info,debug};
use anyhow::{Context, Result};


#[derive(Debug)]
pub struct Heuristic_path{
	pub path: Vec<Indices>,
	pub outpath: Vec<Indices>,
	pub weight: Vec<usize>,
}
impl Heuristic_path{
	pub fn new(path: Vec<Indices>) -> Heuristic_path {
		Heuristic_path {
			path: path,
			outpath: Vec::new(),
			weight: Vec::new(),
		}
	}
}

pub struct Correcter {
	pub name: String,
	pub graph: POAGraph,
}

impl Correcter {
	pub fn new(name: String, graph: POAGraph) -> Correcter {
		Correcter {
			name: name,
			graph: graph,
		}
	}
	// 获取某个节点所有出边的总权重
	fn get_outcome_weight_sum(&self, from_node: NodeIndex<usize>) -> i32 {
		self.graph.edges(from_node).map(|edge| *edge.weight()).sum()
	}
	// 获取从一个节点到另一个节点的边的权重
    fn get_edge_weight(&self, from: NodeIndex<usize>, to: NodeIndex<usize>) -> Result<i32> {
        let edge = self.graph.find_edge(from, to)
            .context(format!("Edge not found between nodes {:?} and {:?}", from, to))?;
        
        let weight = self.graph.edge_weight(edge)
            .context("Edge weight not found")?;
        
        Ok(*weight)
    }
	// 判断是否启用启发式校正
    fn is_heuristic_correction_enabled(&self, from: NodeIndex<usize>, to: NodeIndex<usize>, correct_ratio: f64) -> Result<bool> {
        let from_weight_sum = self.get_outcome_weight_sum(from);
        let edge_weight = self.get_edge_weight(from, to)?;
        // let edge_weight_test = self.get_edge_weight(NodeIndex::new(206), NodeIndex::new(207))?;
		// debug!("(is_heuristic_correction_enabled):edge_weight_test 206->207: {:?}", edge_weight_test);
        let weight_ratio = edge_weight as f64 / from_weight_sum as f64;
		// if weight_ratio == 0.2 {
		// 	use petgraph::Direction;
		// 	for edge in self.graph.edges_directed( NodeIndex::new(206), Direction::Outgoing) {
		// 		let target = edge.target();
		// 		let weight = edge.weight();
		// 		println!(
		// 			"Edge from {:?} to {:?} with weight: {:?}",
		// 			206, target, weight
		// 		);
		// 	}
		// }
        debug!("(is_heuristic_correction_enabled):weight_ratio: {:?}, from: {:?}, from_weight_sum: {:?}, to: {:?}, edge_weight: {:?}", weight_ratio, from, from_weight_sum, to, edge_weight);
        Ok(weight_ratio <= correct_ratio)
    }
	fn heuristic_finder(&self,mut heuristic_path: Heuristic_path, indices: &Indices, correct_fix_ratio: f64) -> (Indices,usize) {
		let mut new_heuristic_path = Vec::new();
		// let mut best_path = Vec::new();
		// let mut min_weight: usize = 0;
		// let mut min_index = 0;
		for i in 0..heuristic_path.path.len() {
			let from_node = heuristic_path.path[i].last().unwrap();
			let outcome_weight_sum = self.get_outcome_weight_sum(*from_node);
			debug!("outcome_weight_sum: {:?}", outcome_weight_sum);
			if outcome_weight_sum == 0 {
				// 到了尾节点，添加进outpath
				heuristic_path.outpath.push(heuristic_path.path[i].clone());
				heuristic_path.weight.push(indices.len());
				continue;
			}
			for edge in self.graph.edges(*from_node) {
				let to_node = edge.target();
				let edge_weight = *edge.weight();
				let weight_ratio = edge_weight as f64 / outcome_weight_sum as f64;
				if weight_ratio >= correct_fix_ratio {
					// 新建一个新的路径，添加到new_heuristic_path
					let mut new_path = heuristic_path.path[i].clone();
					new_path.push(to_node);
					new_heuristic_path.push(new_path.clone());
					// 如果to_node在indices中，添加到new_heuristic_outpath
					if indices.contains(&to_node) {
						heuristic_path.outpath.push(new_path);
						// 找到 to_node 在 indices 中的索引位置
						if let Some(position) = indices.iter().position(|&x| x == to_node) {
							heuristic_path.weight.push(position);
						} else {
						// 如果无法找到 to_node 的位置，根据需求可以做一些处理
						panic!("heuristic_finder to_node not found in indices");
						}
					}
				}
			}
			debug!("(heuristic_finder):heuristic_path: {:?} new_heuristic_outpath: {:?}", heuristic_path, heuristic_path.outpath);
		}

		heuristic_path.path = new_heuristic_path;
		let (best_path, min_weight) = if heuristic_path.outpath.is_empty() {
			// debug!("new_heuristic_outpath is empty");
			// debug!("heuristic_path: {:?}", heuristic_path);
			self.heuristic_finder(heuristic_path, indices, correct_fix_ratio)
		} else {
			let (min_index, &min_weight) = heuristic_path.weight.iter()
				.enumerate()
				.min_by_key(|&(_, weight)| *weight)
				.unwrap(); // Assuming there is at least one element
			let mut best_path: Vec<NodeIndex<usize>> = heuristic_path.outpath[min_index].clone();
			best_path.drain(0..1);
			debug!("heuristic_finder:min_weight: {:?}", min_weight);
			debug!("heuristic_finder:best_path: {:?}", best_path);
			(best_path, min_weight)
		};
		return (best_path,min_weight)
	}


	// 启发式校正
	pub fn heuristic_correct(&self, from: NodeIndex<usize>,indices : &Indices, correct_fix_ratio: f64) -> (Indices,usize) {
		let heuristic_path = Heuristic_path::new(vec![vec![from]]);
		debug!("heuristic_path: {:?} Indices: {:?}", heuristic_path, indices);
		let (best_path,min_index) = self.heuristic_finder(heuristic_path, indices, correct_fix_ratio);
		debug!("best_path: {:?}", best_path);
		return (best_path,min_index);
	}



	// 校正
	pub fn correct_indices(&self, indices: &Indices, qual: &Vec<u8>, correct_ratio: f64, correct_fix_ratio: f64) -> Result<(Indices, Vec<u8>)> {
		let mut correct_indices = Vec::new();
		let mut correct_qual: Vec<u8> = Vec::new();
		let mut i = 1;
		correct_indices.push( indices[0]);
		correct_qual.push(qual[0]);
		debug!("(correct_indices) raw indices(len:{}): {:?}", indices.len(), indices);
		while i < indices.len() {
			let from_node = indices[i-1].clone();
			let to_node = indices[i].clone();
			let heuristic_enabled = self.is_heuristic_correction_enabled(from_node, to_node, correct_ratio)?;
			if heuristic_enabled {
				let (best_path,min_weight) = self.heuristic_correct(from_node, indices, correct_fix_ratio);
				correct_qual.extend(vec![63; best_path.len()]);
				correct_indices.extend(best_path);
				debug!("correct_indices: {:?}", correct_indices);
				debug!("min_weight: {:?}", min_weight);
				// 将循环跳转到i = min_index
				i = min_weight + 1;
			} else {
				correct_indices.push(to_node);
				correct_qual.push(qual[i]);
				i += 1;
			}

		}
		debug!("final correct_indices(len:{}): {:?}", correct_indices.len(), correct_indices);
		Ok((correct_indices, correct_qual))
	}
	// 将indices替换为碱基向量
	pub fn replace_with_node_weights(&self, indices: &Indices) -> Vec<u8> {
		let mut correct_seq: Vec<u8> = Vec::with_capacity(indices.len()); // 预先分配足够的容量
		for &node_index in indices {
			let weight = self.graph.node_weight(node_index)
				.expect("Node index out of bounds or no weight found");
			correct_seq.push(*weight); // 使用 push 方法添加元素
		}
		correct_seq
	}
	// pub fn correct(read) {

	// }

}




#[test]
fn test_main(){
	pretty_env_logger::init();

}

#[test]
fn test_correcter() {
	pretty_env_logger::init();
	debug!("test_correcter");
	let mut graph: Graph<u8, i32, Directed, usize> =
		Graph::with_capacity(10, 9);
	let a = graph.add_node("A".as_bytes()[0]);
	let b = graph.add_node("T".as_bytes()[0]);
	let c = graph.add_node("G".as_bytes()[0]);
	let d = graph.add_node("C".as_bytes()[0]);
	let indices: Vec<NodeIndex<usize>> = vec![a, c, d];
	graph.add_edge(a, b, 20);
	graph.add_edge(a, c, 1);
	graph.add_edge(b, d, 20);
	graph.add_edge(c, d, 1);
	debug!("{:?}", Vec::from(graph.raw_edges()));
	let correcter = Correcter::new("test".to_string(), graph);
	let b_weight_sum = correcter.get_outcome_weight_sum(b);
	debug!("b_weight_sum: {}", b_weight_sum);
	let edge_weight = correcter.get_edge_weight(a, c);
	debug!("edge_weight: {:?}", edge_weight);
	let is_heuristic_correction_enabled = correcter.is_heuristic_correction_enabled(a, c,  0.2);
	debug!("is_heuristic_correction_enabled: {:?}", is_heuristic_correction_enabled);
	// let (correct_indices,correct_qual) = correcter.correct_indices(&indices);
	// debug!("correct_indices: {:?}", correct_indices);
	// let correct_seq = correcter.replace_with_node_weights(&correct_indices);
	// debug!("correct_seq: {:?}", correct_seq);

}
