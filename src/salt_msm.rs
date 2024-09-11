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
        // 确保窗口大小在有效范围内
        assert!(window_size >= 2);
        assert!(window_size < 64);
        // 返回包含窗口大小的新 WnafContext 实例
        Self { window_size }
    }

    /// Generates a precomputed table of multiples of a group element.
    ///
    /// 生成一个群元素倍数的预计算表。
    ///
    /// # Parameters
    /// - `base`: The base group element to generate multiples of.
    ///
    /// # Returns
    /// A vector of group elements representing the precomputed table.
    ///
    /// # 参数
    /// - `base`: 要生成倍数的基群元素。
    ///
    /// # 返回值
    /// 表示预计算表的群元素向量。
    pub fn table<G: Group>(&self, mut base: G) -> Vec<G> {
        // 计算窗口的数量，等价于 ceil(256 / self.window_size)
        let window_count = (256 + self.window_size - 1) / self.window_size;
        // 计算表的大小
        let table_size = window_count * (1 << (self.window_size - 1));
        // 初始化一个具有适当容量的向量
        let mut table = Vec::with_capacity(table_size);
        // 阈值
        let threshold = 1 << (self.window_size - 1);

        // 遍历窗口的数量
        for _ in 0..window_count {
            // 初始化 current_base 为 base
            let mut current_base = base;
            // 遍历窗口中的元素，计算2G,3G,4G,...,2^(window_size-1)G
            //table中保存的为G，2G，3G，4G，...，2^(window_size-1)G
            for _ in 0..threshold {
                // 将 current_base 添加到表中
                table.push(current_base);
                // 更新 current_base，将其加上 base
                current_base += &base;
            }
            // 更新 base 为 current_base，为2^(window_size-1)G+G
            base = current_base;
        }
        // 返回预计算表
        table
    }
    pub fn table2<G: Group>(&self, mut base: G) -> Vec<G> {
        let mut table = Vec::with_capacity(1 << (self.window_size - 1));
        let dbl = base.double();

        for _ in 0..(1 << (self.window_size - 1)) {
            table.push(base);
            base += &dbl;
        }
        table
    }
    /// Computes scalar multiplication of a group element `g` by `scalar`.
    ///
    /// This method uses the wNAF algorithm to perform the scalar
    /// multiplication; first, it uses `Self::table` to calculate an
    /// appropriate table of multiples of `g`, and then uses the wNAF
    /// algorithm to compute the scalar multiple.
    ///
    /// 计算群元素 `g` 与标量 `scalar` 的标量乘法。
    ///
    /// 该方法使用 wNAF 算法执行标量乘法；首先，它使用 `Self::table` 计算 `g` 的适当倍数表，然后使用 wNAF 算法计算标量倍数。
    ///
    /// # Parameters
    /// - `g`: The group element to be multiplied.
    /// - `scalar`: The scalar to multiply the group element by.
    ///
    /// # Returns
    /// The result of the scalar multiplication.
    ///
    /// # 参数
    /// - `g`: 要乘的群元素。
    /// - `scalar`: 要乘以群元素的标量。
    ///
    /// # 返回值
    /// 标量乘法的结果。
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
    ///
    /// 计算群元素与标量 `scalar` 的标量乘法。
    /// `base_table` 包含群元素的预计算倍数；它可以使用 `Self::table` 生成。`scalar` 是 `G::ScalarField` 的元素。
    ///
    /// 如果表太小，则返回 `None`。
    ///
    /// # Parameters
    /// - `base_table`: The precomputed table of group element multiples.
    /// - `scalar`: The scalar to multiply the group element by.
    ///
    /// # Returns
    /// The result of the scalar multiplication, or `None` if the table is too small.
    ///
    /// # 参数
    /// - `base_table`: 群元素倍数的预计算表。
    /// - `scalar`: 要乘以群元素的标量。
    ///
    /// # 返回值
    /// 标量乘法的结果，如果表太小则返回 `None`。
    pub fn mul_with_table1<G: Group>(&self, base_table: &[G], scalar: &G::ScalarField) -> Option<G> {
        // 如果表的大小小于所需的大小，则返回 None
        if 1 << (self.window_size - 1) > base_table.len() {
            return None;
        }

        // 将标量转换为 wNAF 数据
        let data = WnafContext::scalar_to_wnaf_data::<G>(scalar, self.window_size);

        let mut c = 0; // 初始化进位
        let thr = 1 << (self.window_size - 1); // 阈值
        let mut result = G::zero(); // 初始化结果为群的零元素

        // 遍历 wNAF 数据
        for i in 0..data.len() {
            let mut idx = (data[i] + c) as usize; // 计算索引
            if idx == 0 {
                continue; // 如果索引为 0，跳过
            }

            c = 0; // 重置进位
            if idx > thr {
                // 如果索引大于阈值
                idx = (1 << self.window_size) - idx; // 计算新的索引
                if idx != 0 {
                    let a = base_table[idx - 1 + i * thr].clone(); // 获取表中的元素
                    result -= a; // 从结果中减去元素
                }
                c = 1; // 设置进位
            } else {
                let a = base_table[idx - 1 + i * thr]; // 获取表中的元素
                result += a; // 将元素加到结果中
            }
        }
        Some(result) // 返回结果
    }
    pub fn mul_with_table<G: Group>(&self, base_table: &[G], scalar: &G::ScalarField) -> Option<G> {
        // 如果表的大小小于所需的大小，则返回 None
        if 1 << (self.window_size - 1) > base_table.len() {
            return None;
        }

        // 将标量转换为 wNAF 数据
        let wnaf_data = WnafContext::scalar_to_wnaf_data::<G>(scalar, self.window_size);
        let mut carry = 0; // 初始化进位
        let threshold = 1 << (self.window_size - 1); // 阈值
        let mut result = G::zero(); // 初始化结果为群的零元素

        // 遍历 wNAF 数据
        for (i, &wnaf_value) in wnaf_data.iter().enumerate() {
            let mut index = (wnaf_value + carry) as usize; // 计算索引
            if index == 0 {
                continue; // 如果索引为 0，跳过
            }

            carry = 0; // 重置进位
            if index > threshold {
                // 如果索引大于阈值
                index = (1 << self.window_size) - index; // 计算新的索引
                if index != 0 {
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
    #[inline]
    fn scalar_to_wnaf_data<G: Group>(scalar: &G::ScalarField, w: usize) -> Vec<u64> {
        // 将标量转换为 u64 数组
        let source = WnafContext::scalar_to_u64::<G>(scalar);
        // 计算掩码，用于提取 w 位
        let mask = (1 << w) - 1;
        // 初始化存储 wNAF 数据的向量
        let mut data = vec![];
        // 初始化偏移量
        let mut off = w;

        // 遍历源数据
        for i in 0..source.len() {
            // 如果偏移量不等于 w
            let s = if off != w {
                // 计算新的掩码
                let mask = (1 << (w - off)) - 1;
                // 获取上一个数据的索引
                let j = data.len() - 1;
                // 更新上一个数据
                data[j] += (source[i] & mask) << off;
                // 返回新的数据和新的偏移量
                (source[i] >> (w - off), 64 - w + off)
            } else {
                // 否则，直接返回当前数据和偏移量
                (source[i], 64)
            };

            // 按照 w 的步长遍历数据
            for j in (0..s.1).step_by(w) {
                // 提取 w 位数据
                let d = (s.0 >> j) & mask;
                // 将提取的数据添加到 wNAF 数据中
                data.push(d);
                // 更新偏移量
                off = j;
            }
            // 更新偏移量
            off = s.1 - off;
        }

        // 返回 wNAF 数据
        data
    }

    #[inline]
    fn scalar_to_u64<G: Group>(scalar: &G::ScalarField) -> Vec<u64> {
        // 将标量转换为大整数
        let b = scalar.into_bigint();
        // 获取大整数的位数
        let mut num = b.num_bits();
        // 计算需要的 u64 数组长度
        num = if num & 63 == 0 {
            num >> 6 // 如果位数是 64 的倍数，直接除以 64
        } else {
            (num >> 6) + 1 // 否则，除以 64 后加 1
        };

        // 初始化一个具有适当容量的向量
        let mut res = Vec::with_capacity(num as usize);
        // 将大整数的前 num 个 u64 元素添加到向量中
        res.extend_from_slice(&b.as_ref()[..num as usize]);
        // 返回结果向量
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

    /// Generates a vector of random elements based on a given seed.
    ///
    /// 根据给定的种子生成一个随机元素的向量。
    ///
    /// # Parameters
    /// - `num_required_points`: The number of required points.
    /// - `seed`: A seed for the random number generator.
    ///
    /// # Returns
    /// A vector of `Element` instances.
    ///
    /// # 参数
    /// - `num_required_points`: 所需点的数量。
    /// - `seed`: 随机数生成器的种子。
    ///
    /// # 返回值
    /// 一个 `Element` 实例的向量。
    fn generate_random_elements(num_required_points: usize, seed: &'static [u8]) -> Vec<Element> {
        use sha2::{Digest, Sha256};

        let _choose_largest = false;

        // 哈希种子 + i 以获得可能的 x 值
        let hash_to_x = |index: u64| -> Vec<u8> {
            // 创建一个新的 Sha256 哈希器
            let mut hasher = Sha256::new();
            // 更新哈希器的输入数据为种子
            hasher.update(seed);
            // 更新哈希器的输入数据为索引的字节表示
            hasher.update(index.to_be_bytes());
            // 计算哈希值并转换为字节向量
            let bytes: Vec<u8> = hasher.finalize().to_vec();
            // 返回字节向量
            bytes
        };

        // 生成从 0 开始的无限迭代器
        (0u64..)
            // 将迭代器中的每个索引映射为哈希值
            .map(hash_to_x)
            // 尝试将哈希值转换为元素，如果成功则保留
            .filter_map(|hash_bytes| try_reduce_to_element(&hash_bytes))
            // 只保留所需数量的点
            .take(num_required_points)
            // 收集结果为向量
            .collect()
    }
    #[test]
    fn smoke_test_msm() {
        // Generate a base vector of 256 random elements
        let g_base: Vec<_> = generate_random_elements(256, b"eth_verkle_oct_2021").into_iter().collect();

        // Create a precomputed WNAF context with a window size of 10
        let precomp = MSMPrecompWnaf::new(&g_base, 10);

        // Initialize a random number generator with a fixed seed
        let mut rng = ChaCha20Rng::from_seed([0u8; 32]);

        // Create an array of tuples containing random scalars and their indices
        let ax_array: Vec<(Fr, usize)> = (0..32).map(|i| (Fr::rand(&mut rng), i)).collect();

        // Perform MSM using the delta method
        let result = msm_by_delta(&precomp, &ax_array);

        // Assert that the result is not zero
        assert!(!result.is_zero(), "MSM result should not be zero");
    }
    #[test]
    /// Tests the performance of multi-scalar multiplication (MSM) using precomputed WNAF.
    ///
    /// 使用预计算的 WNAF 测试多标量乘法（MSM）的性能。
    fn mertic_msm() {
        // Generate a base vector of 256 random elements
        // 生成一个包含 256 个随机元素的基向量
        let g_base: Vec<_> = generate_random_elements(256, b"eth_verkle_oct_2021").into_iter().collect();

        // Create a precomputed WNAF context with a window size of 10
        // 创建一个窗口大小为 10 的预计算 WNAF 上下文
        let precomp = MSMPrecompWnaf::new(&g_base, 10);

        for num_loads in [32] {
            // Initialize a random number generator with a fixed seed
            // 使用固定种子初始化随机数生成器
            let mut rng = ChaCha20Rng::from_seed([0u8; 32]);

            // Create an array of tuples containing random scalars and their indices
            // 创建一个包含随机标量及其索引的元组数组
            let ax_array: Vec<(Fr, usize)> = (0..num_loads).map(|i| (Fr::rand(&mut rng), i)).collect();

            // Perform MSM using the delta method
            // 使用 delta 方法执行 MSM
            msm_by_delta(&precomp, &ax_array);
        }
    }

    /// Performs multi-scalar multiplication (MSM) using the delta method.
    ///
    /// 使用 delta 方法执行多标量乘法（MSM）。
    ///
    /// # Parameters
    /// - `c`: The precomputed WNAF context.
    /// - `frs`: A slice of tuples containing scalars and their indices.
    ///
    /// # Returns
    /// An `Element` resulting from the MSM.
    ///
    /// # 参数
    /// - `c`: 预计算的 WNAF 上下文。
    /// - `frs`: 包含标量及其索引的元组切片。
    ///
    /// # 返回值
    /// MSM 结果的 `Element`。
    fn msm_by_delta(c: &MSMPrecompWnaf, frs: &[(Fr, usize)]) -> Element {
        let start = Instant::now();
        let e = commit_sparse(c, frs.to_vec());
        let duration = start.elapsed();
        println!("msm_by_delta({}) took: {:?}", frs.len(), duration);
        e
    }

    /// Commits to a sparse set of values using the precomputed WNAF context.
    ///
    /// 使用预计算的 WNAF 上下文提交稀疏值集。
    ///
    /// # Parameters
    /// - `precomp`: The precomputed WNAF context.
    /// - `val_indices`: A vector of tuples containing scalars and their indices.
    ///
    /// # Returns
    /// An `Element` resulting from the commitment.
    ///
    /// # 参数
    /// - `precomp`: 预计算的 WNAF 上下文。
    /// - `val_indices`: 包含标量及其索引的元组向量。
    ///
    /// # 返回值
    /// 提交结果的 `Element`。
    fn commit_sparse(precomp: &MSMPrecompWnaf, val_indices: Vec<(Fr, usize)>) -> Element {
        let mut result = Element::zero();

        for (value, lagrange_index) in val_indices {
            result += precomp.mul_index(value, lagrange_index); // self.scalar_mul(value, lagrange_index)
        }

        result
    }

}