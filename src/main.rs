use crate::avl_tree::SearchTree;
use basic_seed_chainer::avl_tree;
use basic_seed_chainer::seeding_methods;
use basic_seed_chainer::simulation_utils;
use bio::alignment::pairwise::*;
use libwfa::{affine_wavefront::*, bindings::*, mm_allocator::*, penalties::*};
use rayon::prelude::*;
use std::cmp;
use std::str;
use std::sync::Mutex;
use std::time::Instant;

fn main() {
    let print_all_debug = false;
    let use_minimizers = false;
    let use_wfa = true;
    let sketch = false;
    let m_is_substring = false;
    let num_iters = 2000;
    rayon::ThreadPoolBuilder::new()
        .num_threads(20)
        .build_global()
        .unwrap();
    for theta in [0.05] {
        let mut align_cum = vec![];
        let mut recov_cum = vec![];
        let mut chain_cum = vec![];
        let alpha = -((1.0 - theta) as f64).log(4.0);
        let _C = 2. / (1. - 2. * alpha);
        println!("Theta = {}, alpha = {}, C = {}", theta, alpha, _C);
        let ks = 11..17;
        //        let ks = 11..17;
        for k in ks {
            let k = k as usize;
            let d;
            if sketch {
                d = k - 6;
            } else {
                d = 1;
            }
            let w = 10.;
            let window = 2 * d - 1;
            let coeff_scale = 1.0;
            let v = k - d + 1;
            let t = (k - v + 1) / 2 + 1;
            let n = 4.0_f64.powf(k as f64 * (1. - 2. * alpha) / (2. * coeff_scale));
            let n = n as usize;
            let mut s = simulation_utils::gen_rand_string(n);
            let mut s_seeds = seeding_methods::open_sync_seeds(&s, k, v, t).0;
            let zeta_inter = 2. / (1. - 2. * alpha) * 50. / 8.
                * (n as f64).log(4.)
                * (n as f64).ln()
                * ((1. - theta) as f64).powi(-(k as i32));
            let zeta_sketch_inter = zeta_inter * 4. + 2.;
            let zeta_inter = 1. / (6. * zeta_inter);
            //            dbg!(6. * 1. / zeta);
            //            let zeta = 0.00;
            let zeta_sketch_inter = 1. / (6. * zeta_sketch_inter);
            let zeta;
            if sketch {
                zeta = zeta_sketch_inter;
            } else {
                zeta = zeta_inter;
            }
            let zeta = zeta * w;
//            let zeta = 0.;

            let align_times: Mutex<Vec<_>> = Mutex::new(vec![]);
            let chain_times: Mutex<Vec<_>> = Mutex::new(vec![]);
            let recov: Mutex<f64> = Mutex::new(0.);
            let m;
            let start_ind;

            if m_is_substring {
                start_ind = n / 3;
                m = ((n as f64).powf(0.7) + 200.) as usize;
            } else {
                start_ind = 0;
                m = n as usize;
            }
            let offset = 1000000;
            println!(
                "n = {},  m = {}",
                n,
                m
            );
            (0..num_iters)
                .collect::<Vec<usize>>()
                .into_par_iter()
                .for_each(|_l| {
                    let s_repro;
                    if !m_is_substring {
                        s_repro = simulation_utils::gen_rand_string(n);
                    } else {
                        s_repro = simulation_utils::gen_rand_string(2 * k);
                    }
                    let s_p;
                    if m_is_substring {
                        s_p = simulation_utils::gen_mutated_string(
                            &s[start_ind..start_ind + m],
                            theta,
                        )
                        .0;
                    } else {
                        s_p = simulation_utils::gen_mutated_string(
                            &s_repro[start_ind..start_ind + m],
                            theta,
                        )
                        .0;
                    }

                    //                    println!(">str1");
                    //                    println!("{}",str::from_utf8(&s).unwrap());
                    //                    println!(">str2");
                    //                    println!("{}",str::from_utf8(&s_p).unwrap());

                    let now = Instant::now();
                    let s_repro_seeds;
                    let sp_seeds;
                    if !use_minimizers {
                        s_repro_seeds = seeding_methods::open_sync_seeds(&s_repro, k, v, t).0;
                        sp_seeds = seeding_methods::open_sync_seeds(&s_p, k, v, t).0;
                    } else {
                        s_repro_seeds = seeding_methods::minimizer_seeds(&s_repro, window, k).0;
                        sp_seeds = seeding_methods::minimizer_seeds(&s_p, window, k).0;
                    }
                    if print_all_debug {
                        println!("Seeding time {}", now.elapsed().as_secs_f32());
                    }

                    let mut anchors = vec![];
                    anchors.reserve(n);

                    let now = Instant::now();
                    if m_is_substring {
                        for (kmer2, poses2) in sp_seeds.iter() {
                            if s_seeds.contains_key(kmer2) {
                                let poses1 = &s_seeds[kmer2];
                                for pos1 in poses1 {
                                    for pos2 in poses2 {
                                        anchors.push((*pos1, *pos2));
                                    }
                                }
                            }
                        }
                    } else {
                        for (kmer1, poses1) in s_repro_seeds.iter() {
                            if sp_seeds.contains_key(kmer1) {
                                let poses2 = &sp_seeds[kmer1];
                                for pos1 in poses1 {
                                    for pos2 in poses2 {
                                        anchors.push((*pos1, *pos2));
                                    }
                                }
                            }
                        }
                    }
                    anchors.sort();
                    //
                    if print_all_debug {
                        println!("Number of anchors is {}", anchors.len());
                        println!("Anchor finding time {}", now.elapsed().as_secs_f32());
                    }

                    let mut f = vec![w];
                    let mut pointer_array = vec![0; anchors.len()];
                    let mut avl_tree: SearchTree<[usize; 2]> = SearchTree::new();

                    let now = Instant::now();
                    for (i, anchor) in anchors.iter().enumerate() {
                        avl_tree.insert([anchor.1, i]);
                    }
                    avl_tree.update_query_info(
                        [anchors[0].1, 0],
                        w + zeta * (anchors[0].1 + anchors[0].0) as f64,
                        0,
                        100 * anchors[0].0 as usize + offset,
                        100 * anchors[0].1 as usize + offset,
                    );

                    let mut spur = false;
                    for i in 1..anchors.len() {
                        if anchors[i].1 > anchors[i].0 + 1000 || anchors[i].0 > anchors[i].1 + 1000
                        {
                            spur = true;
                        }
                        let best_f_i;
                        let best_j;

                        let (best_score, best_id) = avl_tree.mrq(
                            [0, 0],
                            [anchors[i].1, i],
                            100 * anchors[i].0 as usize + offset,
                            100 * anchors[i].1 as usize + offset,
                        );
                        if spur {
                            //                            dbg!(best_score, &anchors[i], best_id, i);
                        }
                        if best_score == i64::MIN {
                            best_f_i = w;
                            best_j = i;
                        } else {
                            best_j = best_id;
                            best_f_i =
                                best_score as f64 - zeta * (anchors[i].0 + anchors[i].1) as f64 + w;
                        }
                        //                        if best_f_i < 0.0 {
                        //                            best_f_i = 0.0;
                        //                            best_j = i;
                        //                        }
                        f.push(best_f_i);
                        avl_tree.update_query_info(
                            [anchors[i].1, i],
                            best_f_i + zeta * (anchors[i].0 + anchors[i].1) as f64,
                            i,
                            100 * anchors[i].0 as usize + offset,
                            100 * anchors[i].1 as usize + offset,
                        );

                        if best_j != usize::MAX {
                            pointer_array[i] = best_j;
                        }
                    }
                    //            println!("Chaining time {}", now.elapsed().as_secs_f32());
                    chain_times
                        .lock()
                        .unwrap()
                        .push(now.elapsed().as_secs_f32());

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
                    let mut break_length = 0;
                    let mut break_start = false;
                    let mut break_a = 0;
                    let mut break_b = 0;
                    if best_chain.len() > 1 {
                        let range = best_chain[0].1 - best_chain[best_chain.len() - 1].1;
                        for anchors in best_chain.iter() {
                            if anchors.0 != anchors.1 + start_ind {
                                if break_start == false {
                                    break_start = true;
                                    break_a = cmp::min(anchors.0, anchors.1);
                                    break_b = cmp::max(anchors.0, anchors.1);
                                } else {
                                    let min = cmp::min(anchors.0, anchors.1);
                                    let max = cmp::max(anchors.0, anchors.1);
                                    if min < break_a {
                                        break_a = min;
                                    }
                                    if max > break_b {
                                        break_b = max;
                                    }
                                }
                            } else {
                                if break_start {
                                    break_start = false;
                                    break_length += break_b - break_a;
                                }
                            }
                        }
                        let mut rec = recov.lock().unwrap();
                        if break_length > 0 {
                            if print_all_debug{
                                dbg!(break_length);
                            }
                        }
                        if range > break_length {
                            *rec += (range - break_length) as f64;
                        }
                    }

                    let mut gap_intervals = vec![];
                    let mut curr_r_end = 0;
                    let mut curr_q_end = 0;
                    for (iter, anchor) in best_chain.iter().rev().enumerate() {
                        if iter > 0 {
                            let rgap;
                            let qgap;
                            if anchor.0 > curr_r_end + 1 && anchor.1 > curr_q_end + 1 {
                                rgap = (curr_r_end, anchor.0 - 1);
                                qgap = (curr_q_end, anchor.1 - 1);
                                gap_intervals.push((rgap, qgap));
                            } else {
                            }
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
                        if ref_slice.len() > 10000 || query_slice.len() > 10000 {
                            dbg!(gap);
                            continue;
                        }
                        if !use_wfa {
                            let _alignment =
                                aligner.global(&s[gap.0 .0..gap.0 .1], &s_p[gap.1 .0..gap.1 .1]);
                        } else {
                            let mut wavefronts = AffineWavefronts::new_complete(
                                ref_slice.len(),
                                query_slice.len(),
                                &mut penalties,
                                &alloc,
                            );
                            wavefronts.align(ref_slice, query_slice).unwrap();
                        }
                    }

                    align_times
                        .lock()
                        .unwrap()
                        .push(now.elapsed().as_secs_f32());
                });
            let align_mean = align_times.into_inner().unwrap().into_iter().sum::<f32>();
            let chain_mean = chain_times.into_inner().unwrap().into_iter().sum::<f32>();

            let recov_val = *recov.lock().unwrap() / num_iters as f64 / m as f64;
            println!("Value for k {}", k);
            println!("Mean recov for k {} is {}", k, recov_val);
            println!("Mean align {}", align_mean / num_iters as f32);
            println!("Mean chain {}", chain_mean / num_iters as f32);
            align_cum.push(align_mean / num_iters as f32);
            chain_cum.push(chain_mean / num_iters as f32);
            recov_cum.push(recov_val);
        }
        dbg!(align_cum);
        dbg!(chain_cum);
        dbg!(recov_cum);
    }
}
