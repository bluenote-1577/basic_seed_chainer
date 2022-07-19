use crate::avl_tree::SearchTree;
use basic_seed_chainer::avl_tree;
use basic_seed_chainer::seeding_methods;
use basic_seed_chainer::simulation_utils;
use bio::alignment::pairwise::*;
use bio::scores::blosum62;
use block_aligner::cigar::*;
use block_aligner::scan_block::*;
use block_aligner::scores::*;
use libwfa::{affine_wavefront::*, bindings::*, mm_allocator::*, penalties::*};
use std::time::{Duration, Instant};

fn main() {
    let c = 2;
    let ks = 4..14;
    for k in ks {
        let n = 4_usize.pow(k as u32) * c as usize;
        let k = k as usize * c;
        let d = 1;
        let v = k - d + 1;
        let t = (k - v + 1) / 2 + 1;
        let theta = 0.05;
        let s = simulation_utils::gen_rand_string(n);
        let (s_p, pos_mut) = simulation_utils::gen_mutated_string(&s, theta);

        let now = Instant::now();
        let (s_seeds, _num) = seeding_methods::open_sync_seeds(&s, k, v, t);
        let (sp_seeds, _num) = seeding_methods::open_sync_seeds(&s_p, k, v, t);
        println!("Seeding time {}", now.elapsed().as_secs_f32());

        let mut anchors = vec![];

        let now = Instant::now();
        for (kmer1, poses1) in s_seeds {
            if sp_seeds.contains_key(kmer1) {
                let poses2 = &sp_seeds[kmer1];
                for pos1 in poses1 {
                    for pos2 in poses2 {
                        anchors.push((pos1, *pos2));
                    }
                }
            }
        }
        anchors.sort();
        dbg!(anchors.len());
        println!("Anchor finding time {}", now.elapsed().as_secs_f32());

        let mut f = vec![];
        let mut pointer_array = vec![0; anchors.len()];
        let mut avl_tree: SearchTree<[usize; 2]> = SearchTree::new();

        let now = Instant::now();
        for (i, anchor) in anchors.iter().enumerate() {
            avl_tree.insert([anchor.1, i]);
        }
        avl_tree.update_query_info(
            [anchors[0].1, 0],
            1.,
            0,
            anchors[0].0 as usize,
            anchors[0].1 as usize,
        );

        for i in 1..anchors.len() {
            let mut best_f_i = 0. as f64;
            let mut best_j = usize::MAX;
            let mut num_iter = 0;
            let mut max_num_iter = 0;
            let (best_score, best_id) = avl_tree.mrq(
                [0, 0],
                [anchors[i].1, i],
                anchors[i].0 as usize,
                anchors[i].1 as usize,
            );
            if best_score == i64::MIN {
                best_f_i = 0.0;
                best_j = i;
            } else {
                best_j = best_id;
                best_f_i = best_score as f64 + 1.;
                if best_f_i < 0.0 {
                    best_f_i = 0.0;
                    best_j = i;
                }
            }
            f.push(best_f_i);
            avl_tree.update_query_info(
                [anchors[i].1, i],
                best_f_i + 1.,
                i,
                anchors[i].0 as usize,
                anchors[i].1 as usize,
            );

            if best_j != usize::MAX {
                pointer_array[i] = best_j;
            }
        }
        println!("Chaining time {}", now.elapsed().as_secs_f32());

        let mut vec: Vec<_> = f.iter().enumerate().collect();
        vec.sort_by(|(_, v0), (_, v1)| v1.partial_cmp(v0).unwrap());
        let mut curr_i = vec[0].0;
        let mut prev_i = pointer_array[curr_i];
        let mut best_chain = vec![];
        while curr_i != prev_i {
            best_chain.push(anchors[curr_i]);
            curr_i = prev_i;
            prev_i = pointer_array[curr_i];
        }
        best_chain.push(anchors[curr_i]);
        //    dbg!(&best_chain);

        let mut gap_intervals = vec![];
        let mut curr_r_end = 0;
        let mut curr_q_end = 0;
        for anchor in best_chain.iter().rev() {
            let rgap;
            let qgap;
            if anchor.0 > curr_r_end+1 && anchor.1 > curr_q_end+1 {
                rgap = (curr_r_end, anchor.0 - 1);
                qgap = (curr_q_end, anchor.1 - 1);
                gap_intervals.push((rgap, qgap));
            } else {
            }
            curr_r_end = anchor.0 + k - 1;
            curr_q_end = anchor.1 + k - 1;
        }

        let mut penalties = AffinePenalties {
            match_: 0,
            mismatch: 4,
            gap_opening: 6,
            gap_extension: 2,
        };

        //    dbg!(&gap_intervals);
        let score = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
        let now = Instant::now();
        let mut aligner = Aligner::new(-5, -1, &score);
        let alloc = MMAllocator::new(BUFFER_SIZE_8M as u64);
        for gap in gap_intervals {
            let ref_slice = &s[gap.0 .0..gap.0 .1];
            let query_slice = &s[gap.1 .0..gap.1 .1];
//                        let alignment = aligner.global(&s[gap.0 .0..gap.0 .1], &s_p[gap.1 .0..gap.1 .1]);
            let mut wavefronts = AffineWavefronts::new_complete(
                ref_slice.len(),
                query_slice.len(),
                &mut penalties,
                &alloc,
            );
            wavefronts
                .align(ref_slice, query_slice)
                .unwrap();
        }

        //        let block_size = 64;
        //        let gaps = Gaps {
        //            open: -2,
        //            extend: -1,
        //        };
        //        for gap in gap_intervals {
        //            let r = PaddedBytes::from_bytes::<NucMatrix>(&s[gap.0 .0..gap.0 .1], block_size);
        //            let q = PaddedBytes::from_bytes::<NucMatrix>(&s_p[gap.1 .0..gap.1 .1], block_size);
        //
        //            // Align with traceback, but no x drop threshold.
        //            let a = Block::<_, true, false>::align(&q, &r, &NW1, gaps, block_size..=block_size, 0);
        //            let res = a.res();
        //        }

        println!("Alignment time {}", now.elapsed().as_secs_f32());
    }
}
