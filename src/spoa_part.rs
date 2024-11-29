use packed_simd_2::f32x8;
use std::cmp::max;

pub fn custom(&self, query: TextSlice) -> Traceback {
    assert!(self.graph.node_count() != 0);
    let (m, n) = (self.graph.node_count(), query.len());
    let mut max_in_column = vec![(0, 0); n + 1];
    let mut traceback = Traceback::with_capacity(m, n);
    traceback.initialize_scores(self.scoring.gap_open, self.scoring.yclip_prefix);
    let mut topo = Topo::new(&self.graph);

    while let Some(node) = topo.next(&self.graph) {
        let r = self.graph.raw_nodes()[node.index()].weight;
        let r_char = char::from(r);
        let i = node.index() + 1;
        traceback.last = node;
        let prevs: Vec<NodeIndex<usize>> = self.graph.neighbors_directed(node, Incoming).collect();
        traceback.new_row(i, n + 1, self.scoring.gap_open, self.scoring.xclip_prefix, 0, n + 1);

        for query_index in (0..query.len()).step_by(8) {
            let j = query_index + 1;

            // Load query scores into SIMD vector
            let mut query_scores = [0.0; 8];
            for k in 0..8 {
                if query_index + k < query.len() {
                    query_scores[k] = self.scoring.match_fn.score(r, query[query_index + k]) as f32;
                }
            }
            let query_chunk = f32x8::from_slice_unaligned(&query_scores);

            let max_cell = if prevs.is_empty() {
                // No predecessors, only consider match
                let prev_scores = f32x8::splat(traceback.get(0, j - 1).score as f32);
                let scores = prev_scores + query_chunk;

                let mut max_score = scores.extract(0) as i32;
                for k in 1..8 {
                    if query_index + k < query.len() {
                        let score = scores.extract(k) as i32;
                        if score > max_score {
                            max_score = score;
                        }
                    }
                }
                TracebackCell {
                    score: max_score,
                    op: AlignmentOperation::Match(None),
                }
            } else {
                let mut max_cell = TracebackCell {
                    score: MIN_SCORE,
                    op: AlignmentOperation::Match(None),
                };
                for prev_node in &prevs {
                    let i_p = prev_node.index() + 1;
                    let prev_scores_match = f32x8::splat(traceback.get(i_p, j - 1).score as f32);
                    let prev_scores_del = f32x8::splat(traceback.get(i_p, j).score as f32);
                    let scores_match = prev_scores_match + query_chunk;
                    let scores_del = prev_scores_del + f32x8::splat(self.scoring.gap_open as f32);

                    for k in 0..8 {
                        if query_index + k < query.len() {
                            let score_match = scores_match.extract(k) as i32;
                            let score_del = scores_del.extract(k) as i32;
                            let cell_match = TracebackCell {
                                score: score_match,
                                op: AlignmentOperation::Match(Some((i_p - 1, i - 1, r_char))),
                            };
                            let cell_del = TracebackCell {
                                score: score_del,
                                op: AlignmentOperation::Del(Some((i_p - 1, i, r_char))),
                            };
                            max_cell = max(max_cell, max(cell_match, cell_del));
                        }
                    }
                }
                max_cell
            };

            for k in 0..8 {
                if query_index + k < query.len() {
                    let score_ins = traceback.get(i, j - 1 + k).score + self.scoring.gap_open;
                    let cell_ins = TracebackCell {
                        score: score_ins,
                        op: AlignmentOperation::Ins(Some((i - 1, r_char))),
                    };
                    let final_cell = max(max_cell, cell_ins);
                    traceback.set(i, j + k, final_cell);

                    if max_in_column[j + k].0 < final_cell.score {
                        max_in_column[j + k].0 = final_cell.score;
                        max_in_column[j + k].1 = i;
                    }
                }
            }
        }
    }

    let mut max_in_row = (0, 0);
    for (col_index, &(score, col_max_row)) in max_in_column.iter().enumerate() {
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
