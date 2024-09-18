use ark_ec::PrimeGroup;
use ark_ff::{BigInteger, PrimeField, Zero};
use ark_std::vec::Vec;
use crate::Element;

/// A helper type that contains all the context required for computing
/// a window NAF multiplication of a PrimeGroup element by a scalar.
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

    pub fn table<G: PrimeGroup>(&self, mut base: G) -> Vec<G> {
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

    /// Computes scalar multiplication of a PrimeGroup element `g` by `scalar`.
    ///
    /// This method uses the wNAF algorithm to perform the scalar
    /// multiplication; first, it uses `Self::table` to calculate an
    /// appropriate table of multiples of `g`, and then uses the wNAF
    /// algorithm to compute the scalar multiple.
    pub fn mul<G: PrimeGroup>(&self, g: G, scalar: &G::ScalarField) -> G {
        let table = self.table(g);
        self.mul_with_table(&table, scalar).unwrap()
    }

    /// Computes scalar multiplication of a PrimeGroup element by `scalar`.
    /// `base_table` holds precomputed multiples of the PrimeGroup element; it can be
    /// generated using `Self::table`. `scalar` is an element of
    /// `G::ScalarField`.
    ///
    /// Returns `None` if the table is too small.


    pub fn mul_with_table<G: PrimeGroup>(&self, base_table: &[G], scalar: &G::ScalarField) -> Option<G> {
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

    fn scalar_to_wnaf_data<G: PrimeGroup>(scalar: &G::ScalarField, w: usize) -> Vec<i64> {
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
    fn scalar_to_u64<G: PrimeGroup>(scalar: &G::ScalarField) -> Vec<u64> {
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
    pub fn new(t: usize,b: usize) -> Self {
        assert!(b >= 2);
        assert!(b < 64);
        Self { t,b }
    }
    fn fill_window<G: PrimeGroup>(&self, basis: &mut [G], table: &mut [G]) {
        // 如果 basis 数组为空
        if basis.is_empty() {
            // 将 table 数组填充为零元素
            for i in 0..table.len() {
                table[i] = G::zero();
            }
            return;
        }
        let l = table.len();
        // 递归调用 fill_window，传入去掉第一个元素的 basis 数组和 table 数组的前半部分
        self.fill_window(&mut basis[1..], &mut table[.. l/ 2]);
        // 遍历 table 数组的前半部分
        for i in 0..l / 2 {
            // 将 table 数组后半部分的对应元素设置为当前元素与 basis 数组第一个元素的和
            table[l / 2 + i] = table[i] + basis[0];
        }
    }


    //pub fn table<G: PrimeGroup>(&self,mut base: G) -> Vec<G>
    pub fn table<G: PrimeGroup>(&self,base: G)-> Vec<G> {

        let fr_bits =253;
        //b为窗口bit长度，t为预处理的倍数，window_size为每个窗口内的值的所有可能的值
        let window_size = 1 << self.b;
        //用t把标量分为t个部分，如果t=2，那么标量分为2部分，每部分的位数为fr_bits/2，kP=2A1+A0，
        // A1为高位，A0为低位，A1，A0即为划分出来的部分。A1、A0的所有可能的个数为points_per_column，
        // 如t=2，b=5，那么points_per_column=127
        let points_per_column = (fr_bits + self.t - 1) / self.t as usize;
        //窗口的个数
        let window_count = (points_per_column + self.b - 1) / self.b; // 等价于 ceil(256 / self.window_size)
        //let table_size = window_count * (1 << (self.b - 1));
        //let threshold = 1 << (self.b - 1);
        let mut table_basis = Vec::with_capacity(points_per_column);
        let mut dbl = base.clone();
        table_basis.push(base);
        //计算2^0*G,2^2*G,2^4*G,2^6*G,...,2^252*G
        for _i in 1..(points_per_column) {
            for _j in 0..self.t {
                dbl= dbl.double();
            }
            table_basis.push(dbl.clone());
        }
        let mut nn_table = vec![G::zero(); window_count * window_size];
        for i in 0..window_count {
            let w=i;
            let start = w * self.b as usize;
            let mut end = (w + 1) * self.b as usize;
            if end > table_basis.len() {
                end = table_basis.len();
            }
            let window_basis: &mut [G] = &mut table_basis[start..end];
            let mut table = vec![G::zero(); window_size];
            self.fill_window(window_basis, &mut *table);
            for (j, table_element) in table.into_iter().enumerate() {
                nn_table[w * window_size + j] = table_element;
            }
        }
        nn_table
    }

    pub fn mul<G: PrimeGroup>(&self, g: G, scalar: &G::ScalarField) -> G {
        let table = self.table(g);
        self.mul_with_table(&table, scalar).unwrap()
    }
    pub fn mul_with_table<G: PrimeGroup>(&self, base_pre_table: &[G], mon_scalar: &G::ScalarField) -> Option<G> {
        if 1 << (self.b - 1) > base_pre_table.len() {
            return None;
        }

        let window_size = 1 << self.b;
        let mut accum = G::zero();
        let fr_bits = 253;
        let scalar_u64 = WnafGottiContext::scalar_to_u64::<G>(mon_scalar);

        for t_i in 0..self.t {
            if t_i > 0 {
                accum = accum.double();
            }

            let mut curr_window = 0;
            let mut window_scalar = 0;
            let mut window_bit_pos = 0;

            for k in (0..fr_bits).step_by(self.t) {
                let scalar_bit_pos = k + self.t - t_i - 1;
                if scalar_bit_pos < fr_bits && !mon_scalar.is_zero() {
                    let limb = scalar_u64[scalar_bit_pos >> 6];
                    let bit = (limb >> (scalar_bit_pos & 63)) & 1;
                    window_scalar |= (bit as usize) << (self.b - window_bit_pos - 1);
                }

                window_bit_pos += 1;
                if window_bit_pos == self.b {
                    if window_scalar > 0 {
                        let window_precomp = &base_pre_table[curr_window * window_size..(curr_window + 1) * window_size];
                        accum += window_precomp[window_scalar];
                    }
                    curr_window += 1;
                    window_scalar = 0;
                    window_bit_pos = 0;
                }
            }

            if window_scalar > 0 {
                let window_slice = &base_pre_table[curr_window * window_size..(curr_window + 1) * window_size];
                accum += window_slice[window_scalar];
            }
        }

        Some(accum)
    }
    pub fn mul_with_table1<G: PrimeGroup>(&self, base_pre_table: &[G], mon_scalar: &G::ScalarField)-> Option<G> {
        // 检查 base_table 是否太小
        if 1 << (self.b - 1) > base_pre_table.len() {
            //return None;
        }
        // 将标量转换为 wNAF 数据
        //let scalar = mon_scalar.into_bigint();
        //let source = WnafGottiContext::scalar_to_u64::<G>(scalar);
        let window_size = 1 << self.b;
        let mut accum = G::zero();
        let mut add_double_count = 0;
        for t_i in 0..self.t {
            if t_i > 0 {
                accum=accum.double();
                add_double_count += 1;
            }

            let mut curr_window = 0;
            let mut window_scalar = 0;
            let mut window_bit_pos = 0;
            let fr_bits = 253;
            let scalar_u64 = WnafGottiContext::scalar_to_u64::<G>(mon_scalar);
            //let scalar_u64 = crate::salt_msm::WnafContext::scalar_to_wnaf_data::<G>(scalar, self.b);

            let mut k: usize = 0;
            while k < fr_bits {
                let scalar_bit_pos = (k + self.t - t_i - 1) as usize;
                if scalar_bit_pos < fr_bits && !mon_scalar.is_zero() {
                    let limb = scalar_u64[scalar_bit_pos >> 6];
                    let bit = (limb >> (scalar_bit_pos & 63)) & 1;
                    window_scalar |= (bit as usize) << (self.b as usize - window_bit_pos - 1);
                }

                window_bit_pos += 1;
                if window_bit_pos == self.b {
                    if window_scalar > 0 {
                        // let t1 = std::time::Instant::now();
                        // let jiange = t2.elapsed(); // 计算函数执行时间

                        let window_precomp = &base_pre_table[curr_window * window_size..(curr_window + 1) * window_size];
                        //let t2 = std::time::Instant::now();
                        let a = window_precomp[window_scalar];
                        //let takedata = t2.elapsed(); // 计算函数执行时间
                        // let t1 = std::time::Instant::now();
                        // let jiange = t2.elapsed(); // 计算函数执行时间
                        accum+=a;
                        // let add_time = t1.elapsed(); // 计算函数执行时间
                        // println!("add time: {:?}", add_time);
                        //println!("takedata time: {:?}", takedata);
                        // println!("间隔为: {:?}", jiange);
                        // t2 = std::time::Instant::now();
                        add_double_count += 1;
                    }
                    curr_window += 1;

                    window_scalar = 0;
                    window_bit_pos = 0;
                }
                k += self.t;
            }
            if window_scalar > 0 {
                let window_slice = &base_pre_table[curr_window * window_size..(curr_window + 1) * window_size];
                accum+=window_slice[window_scalar];
                add_double_count += 1;
            }

        }


        return Some(accum);
    }
    pub fn msm_with_multiple_tables<G: PrimeGroup>(&self, base_pre_tables: &[&[G]], mon_scalars: &[&G::ScalarField]) -> Option<Vec<G>> {
    if base_pre_tables.len() != mon_scalars.len() {
        return None;
    }

    let mut results = Vec::with_capacity(mon_scalars.len());

    for (base_pre_table, mon_scalar) in base_pre_tables.iter().zip(mon_scalars.iter()) {
        let result = self.mul_with_table1(base_pre_table, mon_scalar)?;
        results.push(result);
    }


    Some(results)
}
    fn scalar_to_u64<G: PrimeGroup>(scalar: &G::ScalarField) -> Vec<u64> {
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

    basic_g.push(Element::prime_subgroup_generator());
    let wnaf_gotti_context = WnafGottiContext::new(2,8);

    // Perform the operation
    let precompute = wnaf_gotti_context.table(basic_g[0].0);
    println!("precompute len: {:?}", precompute.len());

}