use ark_ec::Group; // 引入 ark_ec 库中的 Group trait
use ark_ff::{BigInteger, PrimeField}; // 引入 ark_ff 库中的 BigInteger 和 PrimeField trait
use ark_std::vec::Vec; // 引入 ark_std 库中的 Vec

/// 一个辅助类型，包含计算窗口 NAF 乘法所需的所有上下文
/// 用于将群元素与标量相乘。
pub struct WnafContext {
    pub window_size: usize, // 窗口大小
}

impl WnafContext {
    /// 构造一个新的上下文，窗口大小为 `window_size`。
    ///
    /// # Panic
    ///
    /// 如果 `window_size` 不在 `2 <= window_size < 64` 范围内，则会触发 panic。
    pub fn new(window_size: usize) -> Self {
        assert!(window_size >= 2); // 确保窗口大小大于等于 2
        assert!(window_size < 64); // 确保窗口大小小于 64
        Self { window_size } // 返回新的 WnafContext 实例
    }

    /// 生成一个包含群元素倍数的表。
    pub fn table<G: Group>(&self, mut base: G) -> Vec<G> {
        let window_count = (256 + self.window_size - 1) / self.window_size; // 计算窗口数量，等价于 ceil(256 / self.window_size)
        let table_size = window_count * (1 << (self.window_size - 1)); // 计算表的大小
        let mut table = Vec::with_capacity(table_size); // 预先分配内存
        let threshold = 1 << (self.window_size - 1); // 计算阈值

        for _ in 0..window_count { // 遍历每个窗口
            let mut current_base = base; // 当前基数
            for i in 0..(1 << self.window_size) - 1 { // 遍历窗口内的每个元素
                if i < threshold { // 如果索引小于阈值
                    table.push(current_base.clone()); // 将当前基数添加到表中
                }
                current_base += &base; // 更新当前基数
            }
            base = current_base; // 更新基数
        }

        table // 返回生成的表
    }

    /// 计算群元素 `g` 与标量 `scalar` 的乘积。
    ///
    /// 该方法使用 wNAF 算法进行标量乘法；首先，它使用 `Self::table` 计算
    /// `g` 的适当倍数表，然后使用 wNAF 算法计算标量乘积。
    pub fn mul<G: Group>(&self, g: G, scalar: &G::ScalarField) -> G {
        let table = self.table(g); // 生成倍数表
        self.mul_with_table(&table, scalar).unwrap() // 使用倍数表和标量计算乘积
    }

    /// 计算群元素与标量 `scalar` 的乘积。
    /// `base_table` 包含预计算的群元素倍数；可以使用 `Self::table` 生成。
    /// `scalar` 是 `G::ScalarField` 的一个元素。
    ///
    /// 如果表太小，则返回 `None`。
    pub fn mul_with_table<G: Group>(&self, base_table: &[G], scalar: &G::ScalarField) -> Option<G> {
        if 1 << (self.window_size - 1) > base_table.len() { // 如果表太小
            return None; // 返回 None
        }

        let wnaf_data = WnafContext::scalar_to_wnaf_data::<G>(scalar, self.window_size); // 将标量转换为 wNAF 数据
        let mut carry = 0; // 初始化进位
        let threshold = 1 << (self.window_size - 1); // 计算阈值
        let mut result = G::zero(); // 初始化结果

        for (i, &wnaf_value) in wnaf_data.iter().enumerate() { // 遍历 wNAF 数据
            let mut index = (wnaf_value + carry) as usize; // 计算索引
            if index == 0 { // 如果索引为 0
                continue; // 跳过当前循环
            }

            carry = 0; // 重置进位
            if index > threshold { // 如果索引大于阈值
                index = (1 << self.window_size) - index; // 计算新的索引
                if index != 0 { // 如果索引不为 0
                    let element = base_table[index - 1 + i * threshold].clone(); // 获取表中的元素
                    result -= element; // 从结果中减去元素
                }
                carry = 1; // 设置进位
            } else {
                let element = base_table[index - 1 + i * threshold]; // 获取表中的元素
                result += element; // 将元素加到结果中
            }
        }

        Some(result) // 返回结果
    }

    /// 将标量转换为 wNAF 数据。
    #[inline]
    fn scalar_to_wnaf_data<G: Group>(scalar: &G::ScalarField, w: usize) -> Vec<u64> {
        let source = WnafContext::scalar_to_u64::<G>(scalar); // 将标量转换为 u64 向量
        let mask = (1 << w) - 1; // 计算掩码
        let mut data = vec![]; // 初始化数据向量
        let mut off = w; // 初始化偏移量

        for i in 0..source.len() { // 遍历源数据
            let s = if off != w { // 如果偏移量不等于窗口大小
                let mask = (1 << (w - off)) - 1; // 计算新的掩码
                let j = data.len() - 1; // 获取数据长度
                data[j] += (source[i] & mask) << off; // 更新数据
                (source[i] >> (w - off), 64 - w + off) // 更新源数据和偏移量
            } else {
                (source[i], 64) // 更新源数据和偏移量
            };

            for j in (0..s.1).step_by(w) { // 遍历偏移量
                let d = (s.0 >> j) & mask; // 计算数据
                data.push(d); // 添加数据到向量
                off = j; // 更新偏移量
            }
            off = s.1 - off; // 更新偏移量
        }

        data // 返回数据向量
    }

    /// 将标量转换为 u64 向量。
    #[inline]
    fn scalar_to_u64<G: Group>(scalar: &G::ScalarField) -> Vec<u64> {
        let b = scalar.into_bigint(); // 将标量转换为大整数
        let mut num = b.num_bits(); // 获取大整数的位数
        num = if num & 63 == 0 { // 如果位数是 64 的倍数
            num >> 6 // 右移 6 位
        } else {
            (num >> 6) + 1 // 右移 6 位并加 1
        };

        let mut res = Vec::with_capacity(num as usize); // 预先分配内存
        res.extend_from_slice(&b.as_ref()[..num as usize]); // 将大整数转换为 u64 向量
        res // 返回 u64 向量
    }
}

#[cfg(test)]
mod tests {
    use super::*; // 引入当前模块的所有内容
    use crate::{try_reduce_to_element, msm::MSMPrecompWnaf, Element, Fr}; // 引入 crate 中的相关模块和类型

    use ark_std::rand::SeedableRng; // 引入 ark_std 库中的 SeedableRng trait
    use ark_std::UniformRand; // 引入 ark_std 库中的 UniformRand trait
    use std::time::Instant; // 引入 std 库中的 Instant 类型
    use rand_chacha::ChaCha20Rng; // 引入 rand_chacha 库中的 ChaCha20Rng 类型

    /// 生成随机元素。
    fn generate_random_elements(num_required_points: usize, seed: &'static [u8]) -> Vec<Element> {
        use sha2::{Digest, Sha256}; // 引入 sha2 库中的 Digest 和 Sha256 类型

        let _choose_largest = false; // 初始化选择最大值标志

        // 将种子和索引哈希为可能的 x 值
        let hash_to_x = |index: u64| -> Vec<u8> {
            let mut hasher = Sha256::new(); // 创建新的 Sha256 哈希器
            hasher.update(seed); // 更新哈希器的种子
            hasher.update(index.to_be_bytes()); // 更新哈希器的索引
            let bytes: Vec<u8> = hasher.finalize().to_vec(); // 获取哈希值
            bytes // 返回哈希值
        };

        (0u64..) // 创建一个从 0 开始的无限迭代器
            .map(hash_to_x) // 将迭代器映射为哈希值
            .filter_map(|hash_bytes| try_reduce_to_element(&hash_bytes)) // 过滤并转换为元素
            .take(num_required_points) // 获取所需数量的元素
            .collect() // 收集元素为向量
    }

    #[test]
    fn mertic_msm() {
        let g_base: Vec<_> = generate_random_elements(256, b"eth_verkle_oct_2021").into_iter().collect(); // 生成随机元素

        let precomp = MSMPrecompWnaf::new(&g_base, 10); // 创建新的 MSM 预计算表

        for num_loads in [32] { // 遍历负载数量
            let mut rng = ChaCha20Rng::from_seed([0u8; 32]); // 创建新的随机数生成器
            let ax_array: Vec<(Fr, usize)> = (0..num_loads).map(|i| (Fr::rand(&mut rng), i)).collect(); // 生成随机标量和索引对

            msm_by_delta(&precomp, &ax_array); // 使用预计算表进行 MSM 计算
        }
    }

    /// 使用预计算表进行 MSM 计算。
    fn msm_by_delta(c: &MSMPrecompWnaf, frs: &[(Fr, usize)]) -> Element {
        let start = Instant::now(); // 获取当前时间
        let e = commit_sparse(c, frs.to_vec()); // 提交稀疏数据
        let duration = start.elapsed(); // 计算耗时
        println!("msm_by_delta({}) took: {:?}", frs.len(), duration); // 打印耗时
        e // 返回结果
    }

    /// 稀疏提交。
    fn commit_sparse(precomp: &MSMPrecompWnaf, val_indices: Vec<(Fr, usize)>) -> Element {
        let mut result = Element::zero(); // 初始化结果

        for (value, lagrange_index) in val_indices { // 遍历值和索引对
            result += precomp.mul_index(value, lagrange_index); // 计算并累加结果
        }

        result // 返回结果
    }
}