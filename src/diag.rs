use crate::consts::*;
use crate::model::System;
use nalgebra::{Complex, Const, Matrix2, Matrix6, SymmetricEigen, Vector2, Vector6};

#[derive(Clone)]
pub struct SEud6{
    pub u : SymmetricEigen<Complex<f64>,Const<6>>,
    pub d : SymmetricEigen<Complex<f64>,Const<6>>
}

impl SEud6{
    fn new(u :SymmetricEigen<Complex<f64>,Const<6>>, d : SymmetricEigen<Complex<f64>,Const<6>>,) -> Self{
        SEud6 { u, d }
    }
    pub fn sort(self) -> Self{
        SEud6::new(
            sort_symmetric_eigen_ascending6(self.u), 
            sort_symmetric_eigen_ascending6(self.d),
        )
    }
    pub fn get_cont(op_seud : &Option<Self>) -> &Self{
        match op_seud{
            Some(seud) => seud,
            None => panic!("no seud")
        }
    }
}

#[derive(Clone)]
pub struct SEud2{
    pub u : SymmetricEigen<Complex<f64>,Const<2>>,
    pub d : SymmetricEigen<Complex<f64>,Const<2>>
}

impl SEud2{
    fn new(u :SymmetricEigen<Complex<f64>,Const<2>>, d : SymmetricEigen<Complex<f64>,Const<2>>,) -> Self{
        SEud2 { u, d }
    }
    pub fn sort(self) -> Self{
        SEud2::new(
            sort_symmetric_eigen_ascending2(self.u), 
            sort_symmetric_eigen_ascending2(self.d),
        )
    }
    pub fn get_cont(op_seud : &Option<Self>) -> &Self{
        match op_seud{
            Some(seud) => seud,
            None => panic!("no seud")
        }
    }
}

pub fn diag_66(kk : Vector2<f64>, system : &System) -> SEud6{

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

            let unsorted = SEud6::new(
                SymmetricEigen::new(hamiltonian_u),SymmetricEigen::new(hamiltonian_d)
            );

            unsorted.sort()
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

            let unsorted = SEud6::new(
                SymmetricEigen::new(hamiltonian_u),SymmetricEigen::new(hamiltonian_d)
            );

            unsorted.sort()
        }
        System::UuudddTmdUM(param) => {
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

            let a = 0.01;

            let hamiltonian_u = Matrix6::new(
                mm-j,ed1p,-plu*a,ed3p, mnu*a,ed2p,
                ed1m,mm-j,ed2m, plu,ed3m, -mnu,
                mnu*a,ed2p,mm-j,ed1p,-plu*a, ed3p,
                ed3m, -mnu,ed1m,mm+j,ed2m, plu,
                -plu*a,ed3p, mnu*a,ed2p,mm+j,ed1p,
                ed2m, plu,ed3m,-mnu,ed1m,mm+j
            );

            let hamiltonian_d = Matrix6::new(
                mm+j,ed1p, plu*a,ed3p,-mnu*a,ed2p,
                ed1m,mm+j,ed2m,-plu,ed3m, mnu,
                -mnu*a,ed2p,mm+j,ed1p, plu*a, ed3p,
                ed3m, mnu,ed1m,mm-j,ed2m,-plu,
                plu*a,ed3p,-mnu*a,ed2p,mm-j,ed1p,
                ed2m,-plu,ed3m, mnu,ed1m,mm-j
            );

            let unsorted = SEud6::new(
                SymmetricEigen::new(hamiltonian_u),SymmetricEigen::new(hamiltonian_d)
            );

            unsorted.sort()
        }
        _ => {panic!("invalid system!")}
    }

    

}

pub fn diag_22(kk : Vector2<f64>, system : &System) -> SEud2{

    match system{
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

            let unsorted = SEud2::new(
                SymmetricEigen::new(hamiltonian_u),SymmetricEigen::new(hamiltonian_d)
            );
        

            unsorted.sort()
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

            let unsorted = SEud2::new(
                SymmetricEigen::new(hamiltonian_u),SymmetricEigen::new(hamiltonian_d)
            );
        

            unsorted.sort()
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

            let unsorted = SEud2::new(
                SymmetricEigen::new(hamiltonian_u),SymmetricEigen::new(hamiltonian_d)
            );
        

            unsorted.sort()
        }
        _ => {panic!("invalid system!")}
    }

    

}

//SymmetricEigen<Complex<f64>,Const<6>>を固有値の大きさの昇順に並び替える
pub fn sort_symmetric_eigen_ascending6(
    eigen: SymmetricEigen<Complex<f64>, Const<6>>,
) -> SymmetricEigen<Complex<f64>, Const<6>> {

    let mut indexed_eigenvalues: Vec<(usize, f64)> = eigen
        .eigenvalues
        .iter()
        .enumerate()
        .map(|(i, &val)| (i, val))
        .collect();

    indexed_eigenvalues.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());

    let mut sorted_eigenvalues = Vector6::<f64>::zeros();
    let mut sorted_eigenvectors = Matrix6::<Complex<f64>>::zeros();

    for (new_index, (old_index, _)) in indexed_eigenvalues.iter().enumerate() {
        sorted_eigenvalues[new_index] = eigen.eigenvalues[*old_index];
        sorted_eigenvectors.set_column(new_index, &eigen.eigenvectors.column(*old_index));
    }

    SymmetricEigen {
        eigenvalues: sorted_eigenvalues,
        eigenvectors: sorted_eigenvectors,
    }
}


//SymmetricEigen<Complex<f64>,Const<2>>を固有値の大きさの昇順に並び替える
pub fn sort_symmetric_eigen_ascending2(
    eigen: SymmetricEigen<Complex<f64>, Const<2>>,
) -> SymmetricEigen<Complex<f64>, Const<2>> {

    let mut indexed_eigenvalues: Vec<(usize, f64)> = eigen
        .eigenvalues
        .iter()
        .enumerate()
        .map(|(i, &val)| (i, val))
        .collect();

    indexed_eigenvalues.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());

    let mut sorted_eigenvalues = Vector2::<f64>::zeros();
    let mut sorted_eigenvectors = Matrix2::<Complex<f64>>::zeros();

    for (new_index, (old_index, _)) in indexed_eigenvalues.iter().enumerate() {
        sorted_eigenvalues[new_index] = eigen.eigenvalues[*old_index];
        sorted_eigenvectors.set_column(new_index, &eigen.eigenvectors.column(*old_index));
    }

    SymmetricEigen {
        eigenvalues: sorted_eigenvalues,
        eigenvectors: sorted_eigenvectors,
    }
}