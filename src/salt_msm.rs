use ark_ec::Group;
use ark_ff::{BigInteger, PrimeField};
use ark_std::vec::Vec;

/// A helper type that contains all the context required for computing
/// a window NAF multiplication of a group element by a scalar.
pub struct WnafContext {
    pub window_size: usize,
}

impl WnafContext {
    /// Constructs a new context for a window of size `window_size`.
    ///
    /// # Panics
    ///
    /// This function will panic if not `2 <= window_size < 64`
    pub fn new(window_size: usize) -> Self {
        assert!(window_size >= 2);
        assert!(window_size < 64);
        Self { window_size }
    }

    pub fn table<G: Group>(&self, mut base: G) -> Vec<G> {
        let win_num = if 256 % self.window_size == 0 {
            256 / self.window_size
        } else {
            256 / self.window_size + 1
        };
        let mut table = Vec::with_capacity(win_num * (1 << (self.window_size - 1)));
        let thr = 1 << (self.window_size - 1);

        for _ in 0..win_num {
            let dbl = base;
            for i in 0..(1 << self.window_size) - 1 {
                if i < thr {
                    table.push(base);
                }
                base += &dbl;
            }
        }

        table
    }

    /// Computes scalar multiplication of a group element `g` by `scalar`.
    ///
    /// This method uses the wNAF algorithm to perform the scalar
    /// multiplication; first, it uses `Self::table` to calculate an
    /// appropriate table of multiples of `g`, and then uses the wNAF
    /// algorithm to compute the scalar multiple.
    pub fn mul<G: Group>(&self, g: G, scalar: &G::ScalarField) -> G {
        let table = self.table(g);
        self.mul_with_table(&table, scalar).unwrap()
    }

    /// Computes scalar multiplication of a group element by `scalar`.
    /// `base_table` holds precomputed multiples of the group element; it can be
    /// generated using `Self::table`. `scalar` is an element of
    /// `G::ScalarField`.
    ///
    /// Returns `None` if the table is too small.

    pub fn mul_with_table<G: Group>(&self, base_table: &[G], scalar: &G::ScalarField) -> Option<G> {
        if 1 << (self.window_size - 1) > base_table.len() {
            return None;
        }

        /*let mut result = G::zero();

        for i in 0..256 {
            let a = base_table[i+i*2].clone();
            let start = std::time::Instant::now();
            result += a;
            let d = start.elapsed();
            println!("result took: (+ {:?} ns)", d.as_nanos());
        }*/
        // The modification principle comes from go-ipa(GitHub - crate-crypto/go-ipa: A Go implementation of cryptographic primitives for Verkle Trees)
        //tt_record!("create index start");
        let data = WnafContext::scalar_to_wnaf_data::<G>(scalar, self.window_size);

        let mut c = 0;
        let thr = 1 << (self.window_size - 1);
        let mut result = G::zero();
        //tt_record!("create index end");

        //tt_record!("ecmul start");
        for i in 0..data.len() {
            let mut idx = (data[i] + c) as usize;
            if idx == 0 {
                continue;
            }

            c = 0;
            if idx > thr {
                idx = (1 << self.window_size) - idx;
                if idx != 0 {
                    //tt_record!("ecsub start");
                    //let start = std::time::Instant::now();
                    let a = base_table[idx - 1 + i * thr].clone();
                    //let d = start.elapsed();
                    //println!("precompute get took: (+ {:?} ns)", d.as_nanos());
                    //let start = std::time::Instant::now();
                    result -= a;
                    //let d = start.elapsed();
                    //println!("ecsub took: (+ {:?} ns)", d.as_nanos());
                    //tt_record!("ecsub end");
                }
                c = 1;
            } else {
                //tt_record!("ecadd start");
                let a = base_table[idx - 1 + i * thr];
                //let start = std::time::Instant::now();
                result += a;
                //let d = start.elapsed();
                //println!("ecadd took: (+ {:?} ns)", d.as_nanos());

                //tt_record!("ecadd end");
            }
        }
        //tt_record!("ecmul end");
        Some(result)
    }

    #[inline]
    fn scalar_to_wnaf_data<G: Group>(scalar: &G::ScalarField, w: usize) -> Vec<u64> {
        let source = WnafContext::scalar_to_u64::<G>(scalar);
        let mask = (1 << w) - 1;
        let mut data = vec![];
        let mut off = w;

        for i in 0..source.len() {
            let s = if off != w {
                let mask = (1 << (w - off)) - 1;
                let j = data.len() - 1;
                data[j] += (source[i] & mask) << off;
                (source[i] >> (w - off), 64 - w + off)
            } else {
                (source[i], 64)
            };

            for j in (0..s.1).step_by(w) {
                let d = (s.0 >> j) & mask;
                data.push(d);
                off = j;
            }
            off = s.1 - off;
        }

        data
    }

    #[inline]
    fn scalar_to_u64<G: Group>(scalar: &G::ScalarField) -> Vec<u64> {
        let b = scalar.into_bigint();
        let mut num = b.num_bits();
        num = if num & 63 == 0 {
            num >> 6
        } else {
            (num >> 6) + 1
        };

        let mut res = vec![0u64; num as usize];
        for i in 0..res.len() {
            res[i] = b.as_ref()[i];
        }
        res
    }
}



#[cfg(test)]
mod tests {
    use super::*;
    use crate::{try_reduce_to_element, msm::MSMPrecompWnaf, Element, Fr};

    use ark_std::rand::SeedableRng;
    use ark_std::UniformRand;
    use std::time::Instant;
    use rand_chacha::ChaCha20Rng;

    fn generate_random_elements(num_required_points: usize, seed: &'static [u8]) -> Vec<Element> {
        use sha2::{Digest, Sha256};

        let _choose_largest = false;

        // Hash the seed + i to get a possible x value
        let hash_to_x = |index: u64| -> Vec<u8> {
            let mut hasher = Sha256::new();
            hasher.update(seed);
            hasher.update(index.to_be_bytes());
            let bytes: Vec<u8> = hasher.finalize().to_vec();
            bytes
        };

        (0u64..)
            .map(hash_to_x)
            .filter_map(|hash_bytes| try_reduce_to_element(&hash_bytes))
            .take(num_required_points)
            .collect()
    }

    #[test]
    fn mertic_msm() {
        let g_base: Vec<_> = generate_random_elements(256, b"eth_verkle_oct_2021").into_iter().collect();

        let precomp = MSMPrecompWnaf::new(&g_base, 10);

        for num_loads in [32] {
            //let mut rng = rand::thread_rng();
            let mut rng = ChaCha20Rng::from_seed([0u8; 32]);
            /*let ax_array: Vec<Fr> = (0..256)
                .map(|i| {
                    if i < num_loads {
                        Fr::rand(&mut rng)
                    } else {
                        Fr::from(0)
                    }
                })
                .collect();

            let start = Instant::now();
            msm_by_array(&committer, &ax_array);
            let duration = start.elapsed();
            println!("msm_by_array({}) took: {:?}", num_loads, duration);*/

            let ax_array: Vec<(Fr, usize)> = (0..num_loads).map(|i| (Fr::rand(&mut rng), i)).collect();

            msm_by_delta(&precomp, &ax_array);

        }
    }

    fn msm_by_delta(c: &MSMPrecompWnaf, frs: &[(Fr, usize)]) -> Element {
        let start = Instant::now();
        let e = commit_sparse(c, frs.to_vec());
        let duration = start.elapsed();
        println!("msm_by_delta({}) took: {:?}", frs.len(), duration);
        e
    }

    fn commit_sparse(precomp: &MSMPrecompWnaf, val_indices: Vec<(Fr, usize)>) -> Element {
        let mut result = Element::zero();

        for (value, lagrange_index) in val_indices {
            result += precomp.mul_index(value, lagrange_index) //self.scalar_mul(value, lagrange_index)
        }

        result
    }

}