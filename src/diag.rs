use crate::consts::*;
use crate::model::System;
use nalgebra::{Complex, Const, Matrix2, Matrix6, SymmetricEigen, Vector2, DimMin, Dim, OMatrix, OVector};

#[derive(Clone,Debug)]
pub struct SEud<const N: usize> {
    pub u: SymmetricEigen<Complex<f64>, Const<N>>,
    pub d: SymmetricEigen<Complex<f64>, Const<N>>,
}

impl<const N: usize> SEud<N>
where
    Const<N>: Dim + DimMin<Const<N>, Output = Const<N>>,
      // for matrix
{
    pub fn new(
        u: SymmetricEigen<Complex<f64>, Const<N>>,
        d: SymmetricEigen<Complex<f64>, Const<N>>,
    ) -> Self {
        SEud { u, d }
    }

    pub fn sort(self) -> Self {
        SEud::new(
            sort_symmetric_eigen_ascending(self.u),
            sort_symmetric_eigen_ascending(self.d),
        )
    }
}

#[derive(Clone,Debug)]
pub enum SEudEnum{
    SEud2(SEud<2>),
    SEud6(SEud<6>),
}

impl SEudEnum{
    pub fn sort(self) -> Self{
        match self{
            SEudEnum::SEud2(seud) => {
                SEudEnum::SEud2(seud.sort())
            }
            SEudEnum::SEud6(seud) => {
                SEudEnum::SEud6(seud.sort())
            }
        }
    }
    pub fn is_2(&self) -> &SEud<2>{
        match self {
            SEudEnum::SEud2(seud) => seud,
            _ => panic!("is not 2")
        }
    }
    pub fn is_6(&self) -> &SEud<6>{
        match self {
            SEudEnum::SEud6(seud) => seud,
            _ => panic!("is not62")
        }
    }
    pub fn get_cont(op_seud: &Option<Self>) -> &Self {
        match op_seud {
            Some(seudenum) => seudenum,
            None => panic!("no seud"),
        }
    }
}

pub fn diag(kk : Vector2<f64>, system : &System) -> SEudEnum{
    let ed1p = Complex::exp( I * kk.dot(&D1)) * -T;
    let ed1m = Complex::exp(-I * kk.dot(&D1)) * -T;
    let ed2p = Complex::exp( I * kk.dot(&D2)) * -T;
    let ed2m = Complex::exp(-I * kk.dot(&D2)) * -T;
    let ed3p = Complex::exp( I * kk.dot(&D3)) * -T;
    let ed3m = Complex::exp(-I * kk.dot(&D3)) * -T;

    match system{
        System::Uuuddd(param) => {
            let lambda = param.lambda;
            let j = param.jj;
            let mu = param.mu;

            let plu = {
                Complex::exp( I * kk.dot(&A1)) +
                Complex::exp( I * kk.dot(&A2)) +
                Complex::exp( I * kk.dot(&A3))
            } * I * lambda;
            let mnu = {
                Complex::exp(-I * kk.dot(&A1)) +
                Complex::exp(-I * kk.dot(&A2)) +
                Complex::exp(-I * kk.dot(&A3))
            } * I * lambda;

            let mm = -mu * ONE;

            let hamiltonian_u = Matrix6::new(
                mm-j,ed1p,-plu,ed3p, mnu,ed2p,
                ed1m,mm-j,ed2m,-plu,ed3m, mnu,
                mnu ,ed2p,mm-j,ed1p,-plu, ed3p,
                ed3m, mnu,ed1m,mm+j,ed2m,-plu,
                -plu,ed3p, mnu,ed2p,mm+j,ed1p,
                ed2m,-plu,ed3m, mnu,ed1m,mm+j
            );

            let hamiltonian_d = Matrix6::new(
                mm+j,ed1p, plu,ed3p,-mnu,ed2p,
                ed1m,mm+j,ed2m, plu,ed3m,-mnu,
                -mnu,ed2p,mm+j,ed1p, plu, ed3p,
                ed3m,-mnu,ed1m,mm-j,ed2m, plu,
                plu ,ed3p,-mnu,ed2p,mm-j,ed1p,
                ed2m, plu,ed3m,-mnu,ed1m,mm-j
            );

            let unsorted = SEud::<6>::new(
                SymmetricEigen::new(hamiltonian_u),SymmetricEigen::new(hamiltonian_d)
            );

            SEudEnum::SEud6(unsorted.sort())
        }
        System::UuudddTmd(param) => {
            let lambda = param.lambda;
            let j = param.jj;
            let mu = param.mu;

            let plu = {
                Complex::exp( I * kk.dot(&A1)) +
                Complex::exp( I * kk.dot(&A2)) +
                Complex::exp( I * kk.dot(&A3))
            } * I * lambda;
            let mnu = {
                Complex::exp(-I * kk.dot(&A1)) +
                Complex::exp(-I * kk.dot(&A2)) +
                Complex::exp(-I * kk.dot(&A3))
            } * I * lambda;

            let mm = -mu * ONE;

            let a = 0.0;

            let hamiltonian_u = Matrix6::new(
                mm-j,ed1p,-plu*a,ed3p, mnu*a,ed2p,
                ed1m,mm-j,ed2m,-plu,ed3m, mnu,
                mnu*a,ed2p,mm-j,ed1p,-plu*a, ed3p,
                ed3m, mnu,ed1m,mm+j,ed2m,-plu,
                -plu*a,ed3p, mnu*a,ed2p,mm+j,ed1p,
                ed2m,-plu,ed3m, mnu,ed1m,mm+j
            );

            let hamiltonian_d = Matrix6::new(
                mm+j,ed1p, plu*a,ed3p,-mnu*a,ed2p,
                ed1m,mm+j,ed2m, plu,ed3m,-mnu,
                -mnu*a,ed2p,mm+j,ed1p, plu*a, ed3p,
                ed3m,-mnu,ed1m,mm-j,ed2m, plu,
                plu*a,ed3p,-mnu*a,ed2p,mm-j,ed1p,
                ed2m, plu,ed3m,-mnu,ed1m,mm-j
            );

            let unsorted = SEud::<6>::new(
                SymmetricEigen::new(hamiltonian_u),SymmetricEigen::new(hamiltonian_d)
            );

            SEudEnum::SEud6(unsorted.sort())
        }
        System::Sato(param) => {
            let lambda = param.lambda;
            let j = param.jj;
            let mu = param.mu;

            let diag = {
                2. * lambda * (
                    kk.dot(&A1).sin() +
                    kk.dot(&A2).sin() +
                    kk.dot(&A3).sin() 
                )
            } * ONE;

            let off_diag = {
                Complex::exp( I * kk.dot(&D1)) +
                Complex::exp( I * kk.dot(&D2)) +
                Complex::exp( I * kk.dot(&D3))
            } * -T;

            let mm = -mu * ONE;

            let hamiltonian_u = Matrix2::new(
                diag + mm-j,off_diag,
                off_diag.conj(),diag + mm+j
            );
            let hamiltonian_d = Matrix2::new(
                -diag + mm+j,off_diag,
                off_diag.conj(),-diag + mm-j
            );

            let unsorted = SEud::<2>::new(
                SymmetricEigen::new(hamiltonian_u),SymmetricEigen::new(hamiltonian_d)
            );
        

            SEudEnum::SEud2(unsorted.sort())
        }
        System::Tmd(param) => {
            let lambda = param.lambda;
            let mu = param.mu;

            let diag = {
                2. * lambda * (
                    kk.dot(&A1).sin() +
                    kk.dot(&A2).sin() +
                    kk.dot(&A3).sin() 
                )
            } * ONE;

            let off_diag = {
                Complex::exp( I * kk.dot(&D1)) +
                Complex::exp( I * kk.dot(&D2)) +
                Complex::exp( I * kk.dot(&D3))
            } * -T;

            let mm = -mu * ONE;

            let hamiltonian_u = Matrix2::new(
                diag + mm,off_diag,
                off_diag.conj(), mm
            );
            let hamiltonian_d = Matrix2::new(
                -diag + mm,off_diag,
                off_diag.conj(), mm
            );

            let unsorted = SEud::<2>::new(
                SymmetricEigen::new(hamiltonian_u),SymmetricEigen::new(hamiltonian_d)
            );
        
            SEudEnum::SEud2(unsorted.sort())
        }
        System::SatoTmd(param) => {
            let lambda = param.lambda;
            let j = param.jj;
            let mu = param.mu;

            let diag = {
                2. * lambda * (
                    kk.dot(&A1).sin() +
                    kk.dot(&A2).sin() +
                    kk.dot(&A3).sin() 
                )
            } * ONE;

            let off_diag = {
                Complex::exp( I * kk.dot(&D1)) +
                Complex::exp( I * kk.dot(&D2)) +
                Complex::exp( I * kk.dot(&D3))
            } * -T;

            let mm = -mu * ONE;

            let hamiltonian_u = Matrix2::new(
                diag + mm-j,off_diag,
                off_diag.conj(), mm+j
            );
            let hamiltonian_d = Matrix2::new(
                -diag + mm+j,off_diag,
                off_diag.conj(), mm-j
            );

            let unsorted = SEud::<2>::new(
                SymmetricEigen::new(hamiltonian_u),SymmetricEigen::new(hamiltonian_d)
            );
        
            SEudEnum::SEud2(unsorted.sort())
        }
        _ => {panic!("invalid system!")}
    }
}

pub fn sort_symmetric_eigen_ascending<const N: usize>(
    eigen: SymmetricEigen<Complex<f64>, Const<N>>,
) -> SymmetricEigen<Complex<f64>, Const<N>>
where
    Const<N>: Dim + DimMin<Const<N>, Output = Const<N>>,
{
    let mut indexed_eigenvalues: Vec<(usize, f64)> = eigen
        .eigenvalues
        .iter()
        .enumerate()
        .map(|(i, &val)| (i, val))
        .collect();

    indexed_eigenvalues.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());

    // サイズ N のベクトルと行列を作成
    let mut sorted_eigenvalues = OVector::<f64, Const<N>>::zeros();
    let mut sorted_eigenvectors = OMatrix::<Complex<f64>, Const<N>, Const<N>>::zeros();

    for (new_index, (old_index, _)) in indexed_eigenvalues.iter().enumerate() {
        sorted_eigenvalues[new_index] = eigen.eigenvalues[*old_index];
        sorted_eigenvectors.set_column(new_index, &eigen.eigenvectors.column(*old_index));
    }

    SymmetricEigen {
        eigenvalues: sorted_eigenvalues,
        eigenvectors: sorted_eigenvectors,
    }
}