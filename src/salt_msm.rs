use ark_ec::Group;
use ark_ed_on_bls12_381_bandersnatch::{EdwardsProjective, Fr};
use ark_ff::{BigInteger, PrimeField};
use ark_std::vec::Vec;
use crate::Element;
use crate::msm_gotti::MSMPrecompWnafGotti;

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
        let window_count = (256 + self.window_size - 1) / self.window_size; // 等价于 ceil(256 / self.window_size)
        let table_size = window_count * (1 << (self.window_size - 1));
        let mut table = Vec::with_capacity(table_size);
        let threshold = 1 << (self.window_size - 1);

        for _ in 0..window_count {
            // 初始化 current_base 为 base
            let mut current_base = base;
            // 遍历窗口中的元素,table中为G,2G,3G,...,2^(window_size-1)G
            table.push(current_base);
            for _ in 1..threshold {
                current_base += &base;
                // 将 current_base 添加到表中
                table.push(current_base);
            }
            //2^(window_size-1)G+2^(window_size-1)G=2^(window_size)G
            base = current_base + current_base;
        }
        // 返回预计算表
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


    pub fn mul_with_table_old<G: Group>(&self, base_table: &[G], scalar: &G::ScalarField) -> Option<G> {
        // 检查 base_table 是否太小
        if 1 << (self.window_size - 1) > base_table.len() {
            return None;
        }

        // 将标量转换为 wNAF 数据
        let wnaf_data = WnafContext::scalar_to_wnaf_data_old::<G>(scalar, self.window_size);
        let mut carry = 0;
        let threshold = 1 << (self.window_size - 1);
        let mut result = G::zero();

        // 遍历 wNAF 数据
        for (i, &wnaf_value) in wnaf_data.iter().enumerate() {
            let mut index = (wnaf_value + carry) as usize;
            if index == 0 {
                continue;
            }

            carry = 0;
            if index > threshold {
                // 计算负索引并更新结果
                index = (1 << self.window_size) - index;
                if index != 0 {
                    let element = base_table[index - 1 + i * threshold].clone();
                    result -= element;
                }
                carry = 1;
            } else {
                // 计算正索引并更新结果
                let element = base_table[index - 1 + i * threshold];
                result += element;
            }
        }

        // 返回结果
        Some(result)
    }
    pub fn mul_with_table<G: Group>(&self, base_table: &[G], scalar: &G::ScalarField) -> Option<G> {
        // 检查 base_table 是否太小
        if 1 << (self.window_size - 1) > base_table.len() {
            return None;
        }

        // 将标量转换为 wNAF 数据
        let wnaf_data = WnafContext::scalar_to_wnaf_data::<G>(scalar, self.window_size);

        let pre_comp_size = 1 << (self.window_size - 1);
        let mut result = G::zero();
        // 遍历 wNAF 数据
        for (i, &wnaf_value) in wnaf_data.iter().enumerate() {
            if wnaf_value < 0 {
                let element = base_table[(-wnaf_value - 1) as usize + i * pre_comp_size];
                result -= element;
            } else {
                // 计算正索引并更新结果
                let element = base_table[(wnaf_value - 1) as usize + i * pre_comp_size];
                result += element;
            }
        }
        // 返回结果
        Some(result)
    }

    fn scalar_to_wnaf_data_old<G: Group>(scalar: &G::ScalarField, w: usize) -> Vec<u64> {
        // 将标量转换为 u64 向量，蒙哥马利域转为整数域
        let source = WnafContext::scalar_to_u64::<G>(scalar);
        // mask用来异或取最低的w位
        let mask = (1 << w) - 1;
        // 创建一个空向量来存储 wNAF 数据
        let mut win_data = vec![];
        // 初始化偏移量为窗口大小
        let mut off = w;

        // 遍历源向量
        for i in 0..source.len() {
            // 如果偏移量不等于窗口大小，更新数据向量的最后一个元素
            let s = if off != w {
                let mask = (1 << (w - off)) - 1;
                let j = win_data.len() - 1;
                win_data[j] += (source[i] & mask) << off;
                (source[i] >> (w - off), 64 - w + off)
            } else {
                (source[i], 64)
            };

            // 以窗口大小为步长遍历剩余的位
            for j in (0..s.1).step_by(w) {
                // 提取窗口值并将其推入数据向量
                let d = (s.0 >> j) & mask;
                win_data.push(d);
                off = j;
            }
            // 更新偏移量
            off = s.1 - off;
        }
        win_data

    }
    fn scalar_to_wnaf_data<G: Group>(scalar: &G::ScalarField, w: usize) -> Vec<i64> {
        // 将标量转换为 u64 向量，蒙哥马利域转为整数域
        let source = WnafContext::scalar_to_u64::<G>(scalar);
        // mask用来异或取最低的w位
        let mask = (1 << w) - 1;
        // 创建一个空向量来存储 wNAF 数据
        let mut win_data = vec![];
        // 初始化偏移量为窗口大小
        let mut off = w;

        // 遍历源向量
        for i in 0..source.len() {
            // 如果偏移量不等于窗口大小，更新数据向量的最后一个元素
            let s = if off != w {
                let mask = (1 << (w - off)) - 1;
                let j = win_data.len() - 1;
                win_data[j] += (source[i] & mask) << off;
                (source[i] >> (w - off), 64 - w + off)
            } else {
                (source[i], 64)
            };

            // 以窗口大小为步长遍历剩余的位
            for j in (0..s.1).step_by(w) {
                // 提取窗口值并将其推入数据向量
                let d = (s.0 >> j) & mask;
                win_data.push(d);
                off = j;
            }
            // 更新偏移量
            off = s.1 - off;
        }

        // 把 win_data 变成 <i64>
        let mut data: Vec<i64> = win_data.iter().map(|&x| x as i64).collect();
        let threshold = 1 << (w - 1);

        // 遍历 data，处理进位和负值
        for i in 0..data.len() {
            if data[i] >= threshold {
                data[i] -= 1 << w;
                if i + 1 < data.len() {
                    data[i + 1] += 1;
                } else {
                    data.push(1);
                }
            }
        }

        // 返回 wNAF 数据向量
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

        let mut res = Vec::with_capacity(num as usize);
        res.extend_from_slice(&b.as_ref()[..num as usize]);
        res
    }
}
pub struct WnafGottiContext {
    pub t: usize,
    pub b: usize,
}


impl WnafGottiContext {
    /// Constructs a new context for a window of size `window_size`.
    ///
    /// # Panics
    ///
    /// This function will panic if not `2 <= window_size < 64`
    pub fn new(t: usize,b: usize) -> Self {
        assert!(b >= 2);
        assert!(b < 64);
        Self { t,b }
    }

    pub fn table<G: Group>(&self,mut base: G) -> Vec<G> {

        let fr_bits =253;
        //b为窗口bit长度，t为预处理的倍数，window_size为每个窗口内的值的所有可能的值
        let window_size = 1 << self.b;
        //用t把标量分为t个部分，如果t=2，那么标量分为2部分，每部分的位数为fr_bits/2，kP=2A1+A0，
        // A1为高位，A0为低位，A1，A0即为划分出来的部分。A1、A0的所有可能的个数为points_per_column，
        // 如t=2，b=5，那么points_per_column=127
        let points_per_column = (fr_bits + self.t - 1) / self.t as usize;
        //窗口的个数
        let window_count = (points_per_column + self.b - 1) / self.b; // 等价于 ceil(256 / self.window_size)
        let table_size = window_count * (1 << (self.b - 1));
        let mut table = Vec::with_capacity(table_size);
        let threshold = 1 << (self.b - 1);
        println!("threshold: {:?}", threshold);
        for _ in 0..window_count {
            // 初始化 current_base 为 base
            let mut current_base = base;
            // 遍历窗口中的元素,table中为G,2G,3G,...,2^(window_size-1)G
            table.push(current_base);
            for _ in 1..threshold {
                current_base += &base;
                current_base += &base; // Skip one element
                // 将 current_base 添加到表中
                table.push(current_base);
            }
            //2^(window_size-1)G+2^(window_size-1)G=2^(window_size)G
            base = current_base + current_base;
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
        // 检查 base_table 是否太小
        if 1 << (self.b - 1) > base_table.len() {
            return None;
        }

        // 将标量转换为 wNAF 数据
        let wnaf_data = crate::salt_msm::WnafContext::scalar_to_wnaf_data::<G>(scalar, self.b);

        let pre_comp_size = 1 << (self.b - 1);
        let mut result = G::zero();
        // 遍历 wNAF 数据
        for (i, &wnaf_value) in wnaf_data.iter().enumerate() {
            if wnaf_value < 0 {
                let element = base_table[(-wnaf_value - 1) as usize + i * pre_comp_size];
                result -= element;
            } else {
                // 计算正索引并更新结果
                let element = base_table[(wnaf_value - 1) as usize + i * pre_comp_size];
                result += element;
            }
        }
        // 返回结果
        Some(result)
    }

    fn scalar_to_wnaf_data<G: Group>(scalar: &G::ScalarField, w: usize) -> Vec<i64> {
        // 将标量转换为 u64 向量，蒙哥马利域转为整数域
        let source = crate::salt_msm::WnafContext::scalar_to_u64::<G>(scalar);
        // mask用来异或取最低的w位
        let mask = (1 << w) - 1;
        // 创建一个空向量来存储 wNAF 数据
        let mut win_data = vec![];
        // 初始化偏移量为窗口大小
        let mut off = w;

        // 遍历源向量
        for i in 0..source.len() {
            // 如果偏移量不等于窗口大小，更新数据向量的最后一个元素
            let s = if off != w {
                let mask = (1 << (w - off)) - 1;
                let j = win_data.len() - 1;
                win_data[j] += (source[i] & mask) << off;
                (source[i] >> (w - off), 64 - w + off)
            } else {
                (source[i], 64)
            };

            // 以窗口大小为步长遍历剩余的位
            for j in (0..s.1).step_by(w) {
                // 提取窗口值并将其推入数据向量
                let d = (s.0 >> j) & mask;
                win_data.push(d);
                off = j;
            }
            // 更新偏移量
            off = s.1 - off;
        }

        // 把 win_data 变成 <i64>
        let mut data: Vec<i64> = win_data.iter().map(|&x| x as i64).collect();
        let threshold = 1 << (w - 1);

        // 遍历 data，处理进位和负值
        for i in 0..data.len() {
            if data[i] >= threshold {
                data[i] -= 1 << w;
                if i + 1 < data.len() {
                    data[i + 1] += 1;
                } else {
                    data.push(1);
                }
            }
        }

        // 返回 wNAF 数据向量
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

        let mut res = Vec::with_capacity(num as usize);
        res.extend_from_slice(&b.as_ref()[..num as usize]);
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
#[test]
fn correctness_precompute_one_table_test() {
    // Create a vector of 256 elements, each being a multiple of the prime subgroup generator
    // 创建一个包含 256 个元素的向量，每个元素都是素数子群生成元的倍数
    let mut basic_g = Vec::with_capacity(1);

    basic_g.push(Element::prime_subgroup_generator() * Fr::from(1));
    let wnaf_gotti_context = WnafGottiContext::new(2,5);

    // Perform the operation
    let precompute = wnaf_gotti_context.table(basic_g[0].0);
    println!("precompute len: {:?}", precompute.len());

}