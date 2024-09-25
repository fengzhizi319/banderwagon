use crate::scalar_multi::WnafContext;
use ark_ed_on_bls12_381_bandersnatch::{EdwardsProjective, Fr};
use ark_ff::Zero;
use rayon::prelude::*;

use crate::Element;
#[derive(Clone, Debug)]
pub struct MSMPrecompWnaf {
    pub window_size: usize,
    pub tables: Vec<Vec<EdwardsProjective>>,
}

impl MSMPrecompWnaf {
    pub fn new(bases: &[Element], window_size: usize) -> MSMPrecompWnaf {
        let wnaf_context = WnafContext::new(window_size);

        // parallel generation precompute tables
        MSMPrecompWnaf {
            tables: bases.par_iter().map(|base|{
                wnaf_context.table(base.0)
            }).collect(),
            window_size,
        }
    }

    pub fn fix_ecc_mul(&self, _a: &[Fr]) -> Element {
        use std::str::FromStr;
        let res: Element = Element::zero();

        let scalars = vec![
            Fr::from_str("13108968793781547619861935127046491459309155893440570251786403306729687672800").unwrap()
        ];

        for _i in 0.._a.len() {
            let wnaf_context = WnafContext::new(self.window_size);
            let _res: EdwardsProjective = scalars
                .iter()
                .zip(self.tables.iter())
                .filter(|(scalar, _)| !scalar.is_zero())
                .map(|(scalar, table)| wnaf_context.mul_with_table(table, scalar).unwrap())
                .sum();
        }

        /*for _i in 0.._a.len() {
            res = self.mul_index(scalars[0], 0);
        }*/
        res
    }

    pub fn fixg_ecc_mul(&self, _a: &[Fr]) -> Element {
        let mut res: Element = Element::zero();

        for _i in 0.._a.len() {
            res = self.mul_index(_a[_i], 0);
        }
        res
    }

    pub fn ecc_mul(&self, a: &[Fr]) -> Element {
        let mut res: Element = Element::zero();
        for i in 0..a.len() {
            res = self.mul_index(a[i], i%256);
        }
        res
    }

    pub fn ecc_add(&self, l: usize) -> EdwardsProjective {
        let mut result = EdwardsProjective::zero();
        let vl = self.tables[0].len();
        for i in 0..l {
            result += &self.tables[i%256][i%vl];
        }
        result
    }

    pub fn mul_index(&self, scalar: Fr, index: usize) -> Element {
        let wnaf_context = WnafContext::new(self.window_size);
        Element(
            wnaf_context
                .mul_with_table(&self.tables[index], &scalar)
                .unwrap(),
        )
    }

    pub fn mul(&self, scalars: &[Fr]) -> Element {
        let wnaf_context = WnafContext::new(self.window_size);
        let result: EdwardsProjective = scalars
            .iter()
            .zip(self.tables.iter())
            .filter(|(scalar, _)| !scalar.is_zero())
            .map(|(scalar, table)| wnaf_context.mul_with_table(table, scalar).unwrap())
            .sum();

        Element(result)
    }
    // TODO: This requires more benchmarking and feedback to see if we should
    // TODO put this behind a config flag
    pub fn mul_par(&self, scalars: &[Fr]) -> Element {
        let wnaf_context = WnafContext::new(self.window_size);
        let result: EdwardsProjective = scalars
            .par_iter()
            .zip(self.tables.par_iter())
            .filter(|(scalar, _)| !scalar.is_zero())
            .map(|(scalar, table)| wnaf_context.mul_with_table(table, scalar).unwrap())
            .sum();

        Element(result)
    }
}

#[cfg(test)]
mod tests {
    use std::str::FromStr;
    use ark_ec::CurveGroup;
    use super::*;
    use crate::{multi_scalar_mul, Element};
    use crate::msm_gotti::MSM2DASMPrecompWnafGotti;
    use crate::scalar_multi::MSMExtendPrecompWnaf;

    #[test]
    fn correctness_smoke_test() {
        let mut crs = Vec::with_capacity(256);
        for i in 0..256 {
            crs.push(Element::prime_subgroup_generator() * Fr::from((i + 1) as u64));
        }

        let mut scalars = vec![];
        for i in 0..256 {
            scalars.push(-Fr::from(i + 1));
        }

        let result = multi_scalar_mul(&crs, &scalars);

        let precomp = MSMPrecompWnaf::new(&crs, 12);
        let got_result = precomp.mul(&scalars);
        let got_par_result = precomp.mul_par(&scalars);

        assert_eq!(result, got_result);
        assert_eq!(result, got_par_result);
    }
    #[test]
    fn correctness_2d_asm_msm_mul_one_test() {
        // Create a vector of 256 elements, each being a multiple of the prime subgroup generator
        // 创建一个包含 256 个元素的向量，每个元素都是素数子群生成元的倍数

        let basis_num = 1;
        let mut basic_crs = Vec::with_capacity(basis_num);
        for i in 0..basis_num {
            basic_crs.push(Element::prime_subgroup_generator() * Fr::from((i + 1) as u64));
        }
        let mut scalars = vec![];
        //q-1
        scalars.push(Fr::from_str("13108968793781547619861935127046491459309155893440570251786403306729687672800").unwrap());


        let precompute=MSMExtendPrecompWnaf::new(&basic_crs, 2,4);
        let mem_byte_size=precompute.tables.len()*precompute.tables[0].len()*2*32;
        println!("precompute_size: {:?}", mem_byte_size);
        use std::time::Instant;
        let start = Instant::now();
        let got_result = precompute.mul(&scalars);
        let duration = start.elapsed();
        println!("Time elapsed in mul is: {:?}", duration);

        let affine_result= got_result.0.into_affine();
        let string_x="33549696307925229982445904590536874618633472405590028303463218160177641247209";
        let string_y="19188667384257783945677642223292697773471335439753913231509108946878080696678";
        let x= affine_result.x.to_string();
        let y= affine_result.y.to_string();
        assert_eq!(string_x, x);
        assert_eq!(string_y, y);
        println!("got_result: {:?}", affine_result);
        //let mut scalars = vec![];

    }
}
