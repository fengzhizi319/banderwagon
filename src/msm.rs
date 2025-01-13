use ark_ec::scalar_mul::wnaf::WnafContext;
use ark_ed_on_bls12_381_bandersnatch::{EdwardsProjective, Fr};
use ark_ff::Zero;
use rayon::prelude::*;

use crate::Element;
#[derive(Clone, Debug)]
pub struct MSMPrecompWnaf {
    window_size: usize,
    tables: Vec<Vec<EdwardsProjective>>,
}

impl MSMPrecompWnaf {
    pub fn new(bases: &[Element], window_size: usize) -> MSMPrecompWnaf {
        let wnaf_context = WnafContext::new(window_size);
        let mut tables = Vec::with_capacity(bases.len());

        for base in bases {
            tables.push(wnaf_context.table(base.0));
        }

        MSMPrecompWnaf {
            tables,
            window_size,
        }
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
    fn correctness_smoke_test_compare() {
        let basis_num = 1;
        let mut basic_crs = Vec::with_capacity(basis_num);
        for i in 0..basis_num {
            basic_crs.push(Element::prime_subgroup_generator() * Fr::from((i + 1) as u64));
        }
        let mut scalars = vec![];
        //q-1
        scalars.push(Fr::from_str("13108968793781547619861935127046491459309155893440570251786403306729687672800").unwrap());

        let precompute=MSMPrecompWnaf::new(&basic_crs, 2);

        let got_result = precompute.mul(&scalars);

        let affine_result= got_result.0.into_affine();
        let string_x="33549696307925229982445904590536874618633472405590028303463218160177641247209";
        let string_y="19188667384257783945677642223292697773471335439753913231509108946878080696678";
        let x= affine_result.x.to_string();
        let y= affine_result.y.to_string();
        assert_eq!(string_x, x);
        assert_eq!(string_y, y);
    }
}
