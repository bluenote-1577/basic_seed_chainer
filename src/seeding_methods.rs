use fxhash::hash;
use fxhash::{FxHashMap, FxHashSet};

fn position_min<T: Ord>(slice: &[T]) -> Option<usize> {
    slice
        .iter()
        .enumerate()
        .max_by(|(_, value0), (_, value1)| value1.cmp(value0))
        .map(|(idx, _)| idx)
}

pub fn minimizer_seeds(
    s: &[u8],
    w: usize,
    k: usize,
) -> (FxHashMap<&[u8], FxHashSet<usize>>, usize) {
    let mut minimizer_seeds: FxHashMap<&[u8], FxHashSet<usize>> = FxHashMap::default();
    let mut hashes: Vec<usize> = Vec::new();
    let mut positions_selected: Vec<usize> = Vec::new();

    //hash all k-mers
    for i in 0..s.len() - k + 1 {
        hashes.push(hash(&s[i..i + k]));
    }

    //look at windows
    let mut prev_pos = std::usize::MAX;
    for i in 0..hashes.len() - w + 1 {
        let min_pos = position_min(&hashes[i..i + w]).unwrap() + i;
        if prev_pos == min_pos {
            continue;
        } else {
            let pos_vec = minimizer_seeds
                .entry(&s[min_pos..min_pos + k])
                .or_insert(FxHashSet::default());
            pos_vec.insert(min_pos);
            positions_selected.push(min_pos);
            prev_pos = min_pos;
        }
    }

    return (minimizer_seeds, positions_selected.len());
}

pub fn open_sync_seeds(
    string: &[u8],
    k: usize,
    s: usize,
    t: usize,
) -> (FxHashMap<&[u8], FxHashSet<usize>>, usize) {
    let mut syncmer_seeds: FxHashMap<&[u8], FxHashSet<usize>> = FxHashMap::default();
    let mut hashes: Vec<usize> = Vec::new();
    let mut positions_selected: Vec<usize> = Vec::new();
    //hash all s-mers
    for i in 0..string.len() - s + 1 {
        hashes.push(hash(&string[i..i + s]));
    }

    let w = k - s + 1;
    for i in 0..hashes.len() - w + 1 {
        let min_pos = position_min(&hashes[i..i + w]).unwrap() + i;
        if min_pos - i == t - 1 {
            let pos_vec = syncmer_seeds
                .entry(&string[i..i + k])
                .or_insert(FxHashSet::default());
            pos_vec.insert(i);
            positions_selected.push(i)
        }
    }

    return (syncmer_seeds, positions_selected.len());
}

fn get_charged_contexts(string: &[u8], w: usize, k: usize) -> FxHashSet<usize> {
    let mut hashes: Vec<usize> = Vec::new();
    let mut charged_contexts = FxHashSet::default();
    for i in 0..string.len() - k + 1 {
        hashes.push(hash(&string[i..i + k]));
    }

    for i in 0..hashes.len() - (w + 1) + 1 {
        let min_pos = position_min(&hashes[i..i + (w + 1)]).unwrap();
        if min_pos == 0 || min_pos == w {
            charged_contexts.insert(i);
        }
    }

    //dbg!(charged_contexts.len());
    return charged_contexts;
}

fn get_syncmer_contexts(string: &[u8], s: usize, k: usize, t: usize) -> FxHashSet<usize> {
    let mut hashes: Vec<usize> = Vec::new();
    let mut charged_contexts = FxHashSet::default();
    for i in 0..string.len() - s + 1 {
        hashes.push(hash(&string[i..i + s]));
    }

    let w = k - s + 1;
    for i in 0..hashes.len() - w + 1 {
        let min_pos = position_min(&hashes[i..i + w]).unwrap();
        if min_pos == t - 1 {
            charged_contexts.insert(i);
        }
    }

    //dbg!(charged_contexts.len());
    return charged_contexts;
}

pub fn custom_words_seeds_5_4(
    string: &[u8],
    n: usize,
    k: usize,
) -> (FxHashMap<&[u8], FxHashSet<usize>>, usize) {
    let mut W = Vec::<[u8; 5]>::new();
    let A = 65;
    let T = 84;
    let C = 67;
    let G = 71;
    let r = 1;
    let y = 2;
    W.push([r, r, r, r, y]);
    W.push([r, r, y, r, r]);
    W.push([r, y, r, y, r]);
    W.push([r, y, y, r, r]);
    W.push([r, y, y, r, y]);
    W.push([r, y, y, y, y]);
    W.push([y, y, y, r, r]);
    W.push([y, y, y, r, y]);
    let mut seeds = FxHashMap::default();
    let mut positions_selected: Vec<usize> = Vec::new();

    let mut cur_length = 0;
    let mut vec_mer = vec![];
    for i in 0..string.len() - k + 1 {
        cur_length += 1;
        if string[i] == A || string[i] == G {
            vec_mer.push(r);
        } else {
            vec_mer.push(y);
        }
        if cur_length < 5 {
            continue;
        }
        let last_word = &vec_mer[i - 4..i + 1];
        let last_word_sl = [
            last_word[0],
            last_word[1],
            last_word[2],
            last_word[3],
            last_word[4],
        ];
        if W.contains(&last_word_sl) {
            positions_selected.push(i - 4);
            let pos_vec = seeds
                .entry(&string[i - 4..i - 4 + k])
                .or_insert(FxHashSet::default());
            pos_vec.insert(i - 4);
        }
    }

    return (seeds, positions_selected.len());
}

pub fn custom_words_seeds_8_8(
    string: &[u8],
    n: usize,
    k: usize,
) -> (FxHashMap<&[u8], FxHashSet<usize>>, usize) {
    let mut W = Vec::<[u8; 8]>::new();
    let A = 65;
    let T = 84;
    let C = 67;
    let G = 71;
    let r = 1;
    let y = 2;
    let mer_size = 8;
    let words_set = vec![
        "rrrrrrry", "rryrrryy", "ryrrrryr", "ryrrrryy", "ryrrryry", "yrrrrrry", "yrrrrryr",
        "yrrrrryy", "yrrryrry", "yrryrryr", "yrryrryy", "yryrrryy", "yryrryry", "yryrryyr",
        "yryrryyy", "yryryryy", "yyrrrryr", "yyrrrryy", "yyrrryry", "yyrrryyr", "yyrrryyy",
        "yyrryryr", "yyrryryy", "yyrryyry", "yyrryyyr", "yyrryyyy", "yyryryyr", "yyryryyy",
        "yyryyryy", "yyyryyyr", "yyyryyyy", "yyyyyyyr",
    ];

    for w in words_set{
        let mut topush = [0,0,0,0,0,0,0,0];
        for (i,c) in w.chars().enumerate(){
            if c == 'r'{
                topush[i] = 1;
            }
            else{
                topush[i] = 2;
            }
        }
        W.push(topush);
    }

    let mut seeds = FxHashMap::default();
    let mut positions_selected: Vec<usize> = Vec::new();

    let mut cur_length = 0;
    let mut vec_mer = vec![];
    for i in 0..string.len() - k + 1 {
        cur_length += 1;
        if string[i] == A || string[i] == G {
            vec_mer.push(r);
        } else {
            vec_mer.push(y);
        }
        if cur_length < mer_size {
            continue;
        }
        let last_word = &vec_mer[i + 1 - mer_size..i + 1];
        let last_word_sl = [
            last_word[0],
            last_word[1],
            last_word[2],
            last_word[3],
            last_word[4],
            last_word[5],
            last_word[6],
            last_word[7],
        ];
        if W.contains(&last_word_sl) {
            positions_selected.push(i + 1 - mer_size);
            let pos_vec = seeds
                .entry(&string[i + 1 - mer_size..i + 1 - mer_size + k])
                .or_insert(FxHashSet::default());
            pos_vec.insert(i + 1 - mer_size);
        }
    }

    return (seeds, positions_selected.len());
}

pub fn custom_words_seeds_6_4(
    string: &[u8],
    n: usize,
    k: usize,
) -> (FxHashMap<&[u8], FxHashSet<usize>>, usize) {
    let mut W = Vec::<[u8; 6]>::new();
    let A = 65;
    let T = 84;
    let C = 67;
    let G = 71;
    let r = 1;
    let y = 2;
    let mer_size = 6;
    //    rrrrry, rryrry, rryryy, ryrrrr,
    //ryrrry, ryryry, ryyrrr, ryyrry,
    //ryyryr, ryyryy, ryyyry, ryyyyy,
    //yryrry, yyyrrr, yyyrry, yyyyry
    W.push([r, r, r, r, r, y]);
    W.push([r, r, y, r, r, y]);
    W.push([r, r, y, r, y, y]);
    W.push([r, y, r, r, r, r]);

    W.push([r, y, r, r, r, y]);
    W.push([r, y, r, y, r, y]);
    W.push([r, y, y, r, r, r]);
    W.push([r, y, y, r, r, y]);

    W.push([r, y, y, r, y, r]);
    W.push([r, y, y, r, y, y]);
    W.push([r, y, y, y, r, y]);
    W.push([r, y, y, y, y, y]);

    W.push([y, r, y, r, r, y]);
    W.push([y, y, y, r, r, r]);
    W.push([y, y, y, r, r, y]);
    W.push([y, y, y, y, r, y]);
    let mut seeds = FxHashMap::default();
    let mut positions_selected: Vec<usize> = Vec::new();

    let mut cur_length = 0;
    let mut vec_mer = vec![];
    for i in 0..string.len() - k + 1 {
        cur_length += 1;
        if string[i] == A || string[i] == G {
            vec_mer.push(r);
        } else {
            vec_mer.push(y);
        }
        if cur_length < mer_size {
            continue;
        }
        let last_word = &vec_mer[i + 1 - mer_size..i + 1];
        let last_word_sl = [
            last_word[0],
            last_word[1],
            last_word[2],
            last_word[3],
            last_word[4],
            last_word[5],
        ];
        if W.contains(&last_word_sl) {
            positions_selected.push(i + 1 - mer_size);
            let pos_vec = seeds
                .entry(&string[i + 1 - mer_size..i + 1 - mer_size + k])
                .or_insert(FxHashSet::default());
            pos_vec.insert(i + 1 - mer_size);
        }
    }

    return (seeds, positions_selected.len());
}

pub fn words_seeds(
    string: &[u8],
    n: usize,
    k: usize,
) -> (FxHashMap<&[u8], FxHashSet<usize>>, usize) {
    let mut seeds = FxHashMap::default();
    let mut positions_selected: Vec<usize> = Vec::new();
    let a = 65;
    for i in 0..string.len() - k + 1 {
        if string[i] == a {
            let mut prefix_good = true;
            for j in 1..n + 1 {
                if string[i + j] == a {
                    prefix_good = false;
                    break;
                }
            }
            if prefix_good {
                positions_selected.push(i);
                let pos_vec = seeds
                    .entry(&string[i..i + k])
                    .or_insert(FxHashSet::default());
                pos_vec.insert(i);
            }
        }
    }

    return (seeds, positions_selected.len());
}

pub fn miniception_seeds(
    string: &[u8],
    w: usize,
    k: usize,
    k_0: usize,
) -> (FxHashMap<&[u8], FxHashSet<usize>>, usize) {
    let charged_contexts = get_charged_contexts(string, k - k_0, k_0);
    //let charged_contexts = get_syncmer_contexts(string, k_0, k, (k-k_0+1)/2 + 1 as usize);
    let mut minimizer_seeds: FxHashMap<&[u8], FxHashSet<usize>> = FxHashMap::default();
    let mut hashes: Vec<usize> = Vec::new();
    let mut positions_selected: Vec<usize> = Vec::new();

    //hash all k-mers
    for i in 0..string.len() - k + 1 {
        hashes.push(hash(&string[i..i + k]));
    }

    //look at windows
    let mut prev_pos = std::usize::MAX;
    for i in 0..hashes.len() - w + 1 {
        let mut context_hashes = Vec::new();
        let mut non_context_hashes = Vec::new();
        for j in i..(i + w) {
            if charged_contexts.contains(&j) {
                context_hashes.push((hashes[j], j));
            } else {
                non_context_hashes.push((hashes[j], j));
            }
        }

        let hash_set;
        if context_hashes.len() != 0 {
            hash_set = context_hashes;
        } else {
            hash_set = non_context_hashes;
        }

        let min_pos = position_min(&hash_set).unwrap();
        let min_pos = hash_set[min_pos].1;
        if prev_pos == min_pos {
            continue;
        } else {
            let pos_vec = minimizer_seeds
                .entry(&string[min_pos..min_pos + k])
                .or_insert(FxHashSet::default());
            pos_vec.insert(min_pos);
            positions_selected.push(min_pos);
            prev_pos = min_pos;
        }
    }

    return (minimizer_seeds, positions_selected.len());
}
