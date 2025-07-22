use crate::consts::*;
use crate::diag::{diag_uuuddd2,SEud};
use crate::dv2::DV2;

use nalgebra::{Complex, Vector2, Vector6};
use std::fs;
use std::fmt::Write;

pub fn calculate_berry_curvature_z_fukui(export_dat : bool, graph_mesh : usize, lambda : f64, jj : f64, mu : f64, comment : &str) -> Vec<Vec<f64>>{
    
    let mesh_kx = graph_mesh;
    let mesh_ky  = graph_mesh;

    let dv2_k1 = DV2::from_car(KPPKS) - DV2::from_car(-KINKS);
    let dv2_k2 = DV2::from_car(KP_KS) - DV2::from_car(-KINKS);

    let mut seud_mat = vec![vec![None ; mesh_ky + 1] ; mesh_kx + 1];
    let mut berry_curvature_z = vec![vec![0f64 ; mesh_ky] ; mesh_kx];


    for i in 0..(mesh_kx+1){
        for j in 0..(mesh_ky+1){
            let if64 = i as f64 / mesh_kx as f64;
            let jf64 = j as f64 / mesh_ky as f64;

            let kk_dv2 = dv2_k1 * if64 + dv2_k2 * jf64 + DV2::from_car(-KINKS);

            let seud = diag_uuuddd2(kk_dv2.to_car(),lambda,jj,mu).sort();

            seud_mat[i][j] = Some(seud);
        }
    }

    let mut chern = 0.;
    let mut file_str = format!("kx ky berry \ngraphmesh {} lambda {} jj {} mu {}",graph_mesh,lambda, jj, mu);

    for i in 0..mesh_kx{
        for j in 0..mesh_ky{
            let if64 = i as f64 / mesh_kx as f64;
            let jf64 = j as f64 / mesh_ky as f64;

            let mut kk_dv2 = dv2_k1 * if64 + dv2_k2 * jf64 + DV2::from_car(-KINKS);

            if point_in_triangle_simple(kk_dv2.to_car(), KINKS, 2. * KINKS, KPPKS){
                kk_dv2 = kk_dv2 - dv2_k1;
            }
            else if point_in_triangle_simple(kk_dv2.to_car(), KINKS, 2. * KINKS, KP_KS){
                kk_dv2 = kk_dv2 - dv2_k2;
            }
            

            let kx = kk_dv2.to_car().x;
            let ky = kk_dv2.to_car().y;

            let seud = SEud::get_cont(seud_mat[i][j]);

            let se_u = seud.u;
            let se_d = seud.d;

            let mut berry_sum = 0.;

            for eigen_i in 0..se_u.eigenvalues.len(){
                if se_u.eigenvalues[eigen_i] < mu{
                    let u00 : Vector6<Complex<f64>> = se_u.eigenvectors.column(eigen_i).into();
                    let u10 : Vector6<Complex<f64>> = SEud::get_cont(seud_mat[i + 1][j]).u.eigenvectors.column(eigen_i).into();
                    let u11 : Vector6<Complex<f64>> = SEud::get_cont(seud_mat[i + 1][j + 1]).u.eigenvectors.column(eigen_i).into();
                    let u01 : Vector6<Complex<f64>> = SEud::get_cont(seud_mat[i][j + 1]).u.eigenvectors.column(eigen_i).into();

                    let berry = {
                        let u1 = u00.dotc(&u10);
                        let u2 = u10.dotc(&u11);
                        let u3 = u11.dotc(&u01);
                        let u4 = u01.dotc(&u00);

                        (u1 * u2 * u3 * u4).arg()
                    };

                    berry_sum += berry;
                }
            }

            for eigen_i in 0..se_d.eigenvalues.len(){
                if se_d.eigenvalues[eigen_i] < mu{
                    let u00 : Vector6<Complex<f64>> = se_d.eigenvectors.column(eigen_i).into();
                    let u10 : Vector6<Complex<f64>> = SEud::get_cont(seud_mat[i + 1][j]).d.eigenvectors.column(eigen_i).into();
                    let u11 : Vector6<Complex<f64>> = SEud::get_cont(seud_mat[i + 1][j + 1]).d.eigenvectors.column(eigen_i).into();
                    let u01 : Vector6<Complex<f64>> = SEud::get_cont(seud_mat[i][j + 1]).d.eigenvectors.column(eigen_i).into();

                    let berry = {
                        let u1 = u00.dotc(&u10);
                        let u2 = u10.dotc(&u11);
                        let u3 = u11.dotc(&u01);
                        let u4 = u01.dotc(&u00);

                        (u1 * u2 * u3 * u4).arg()
                    };

                    berry_sum += berry;
                }
            }

            if export_dat{
                write!(file_str, "\n{} {} {}", kx, ky, berry_sum).unwrap();
            }

            chern += berry_sum;

            berry_curvature_z[i][j] = berry_sum;
        

        }
    }

    if export_dat{
        fs::write(
            format!(
                "./output/berry/dats/berry_fukui_mesh_{}_lambda_{}_j{}_mu{}_{}.dat"
                ,graph_mesh
                ,lambda.to_string().replace('.', "p")
                ,jj.to_string().replace('.', "p")
                ,mu.to_string().replace('.', "p")
                ,comment
            ), 
            file_str
        ).unwrap();
    }

    println!("{}",chern / 2. / PI);

    berry_curvature_z

}

#[derive(Debug, Clone, Copy)]
struct BerryEigen{
    berry : f64,
    eigen : f64,
    spin : Spin,
}

impl BerryEigen{
    fn ini() -> Self{
        BerryEigen { berry: 0.0, eigen: 0.0, spin : Spin::Undefined }
    }
    fn new(berry : f64, eigen : f64, spin : Spin) -> Self{
        BerryEigen { berry, eigen , spin}
    }
}

#[derive(Debug, Clone, Copy)]
pub enum Spin{
    U,
    D,
    Undefined,
}

impl Spin{
    fn to_u8(self) -> u8{
        match self{
            Self::U => {1}
            Self::D => {0}
            Self::Undefined => {100}
        }
    }
}

pub fn calculate_band_info_all_band(export_dat : bool, graph_mesh : usize, lambda : f64, jj : f64, mu : f64, comment : &str){
    
    let mesh_kx = graph_mesh;
    let mesh_ky  = graph_mesh;

    let dv2_k1 = DV2::from_car(KPPKS) - DV2::from_car(-KINKS);
    let dv2_k2 = DV2::from_car(KP_KS) - DV2::from_car(-KINKS);

    let mut seud_mat = vec![vec![None ; mesh_ky + 1] ; mesh_kx + 1];

    for i in 0..(mesh_kx+1){
        for j in 0..(mesh_ky+1){
            let if64 = i as f64 / mesh_kx as f64;
            let jf64 = j as f64 / mesh_ky as f64;

            let kk_dv2 = dv2_k1 * if64 + dv2_k2 * jf64 + DV2::from_car(-KINKS);

            let seud = diag_uuuddd2(kk_dv2.to_car(),lambda,jj,mu).sort();

            seud_mat[i][j] = Some(seud);
        }
    }

    let mut file_str = format!("kx ky b1 b2 b3 b4 b5 b6 e1 e2 e3 e4 e5 s1 s2 s3 s4 s5 s6\ngraphmesh {} lambda {} jj {} mu {}",graph_mesh,lambda, jj, mu);

    for i in 0..mesh_kx{
        for j in 0..mesh_ky{
            let if64 = i as f64 / mesh_kx as f64;
            let jf64 = j as f64 / mesh_ky as f64;

            let mut kk_dv2 = dv2_k1 * if64 + dv2_k2 * jf64 + DV2::from_car(-KINKS);

            if point_in_triangle_simple(kk_dv2.to_car(), KINKS, 2. * KINKS, KPPKS){
                kk_dv2 = kk_dv2 - dv2_k1;
            }
            else if point_in_triangle_simple(kk_dv2.to_car(), KINKS, 2. * KINKS, KP_KS){
                kk_dv2 = kk_dv2 - dv2_k2;
            }
            

            let kx = kk_dv2.to_car().x;
            let ky = kk_dv2.to_car().y;

            let seud = SEud::get_cont(seud_mat[i][j]);

            let se_u = seud.u;
            let se_d = seud.d;

            let mut berry_s = [BerryEigen::ini();12];

            for eigen_i in 0..se_u.eigenvalues.len(){
                {
                    let u00 : Vector6<Complex<f64>> = se_u.eigenvectors.column(eigen_i).into();
                    let u10 : Vector6<Complex<f64>> = SEud::get_cont(seud_mat[i + 1][j]).u.eigenvectors.column(eigen_i).into();
                    let u11 : Vector6<Complex<f64>> = SEud::get_cont(seud_mat[i + 1][j + 1]).u.eigenvectors.column(eigen_i).into();
                    let u01 : Vector6<Complex<f64>> = SEud::get_cont(seud_mat[i][j + 1]).u.eigenvectors.column(eigen_i).into();

                    let berry = {
                        let u1 = u00.dotc(&u10);
                        let u2 = u10.dotc(&u11);
                        let u3 = u11.dotc(&u01);
                        let u4 = u01.dotc(&u00);

                        (u1 * u2 * u3 * u4).arg()
                    };

                    berry_s[eigen_i] = BerryEigen::new(berry,se_u.eigenvalues[eigen_i],Spin::U);
                }
            }

            for eigen_i in 0..se_d.eigenvalues.len(){
                {
                    let u00 : Vector6<Complex<f64>> = se_d.eigenvectors.column(eigen_i).into();
                    let u10 : Vector6<Complex<f64>> = SEud::get_cont(seud_mat[i + 1][j]).d.eigenvectors.column(eigen_i).into();
                    let u11 : Vector6<Complex<f64>> = SEud::get_cont(seud_mat[i + 1][j + 1]).d.eigenvectors.column(eigen_i).into();
                    let u01 : Vector6<Complex<f64>> = SEud::get_cont(seud_mat[i][j + 1]).d.eigenvectors.column(eigen_i).into();

                    let berry = {
                        let u1 = u00.dotc(&u10);
                        let u2 = u10.dotc(&u11);
                        let u3 = u11.dotc(&u01);
                        let u4 = u01.dotc(&u00);

                        (u1 * u2 * u3 * u4).arg()
                    };

                    berry_s[eigen_i + 6] = BerryEigen::new(berry,se_d.eigenvalues[eigen_i],Spin::D);
                }
            }

           berry_s.sort_by(|a, b| a.eigen.partial_cmp(&b.eigen).unwrap());


            if export_dat{
                write!(file_str, "\n{} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {}", 
                    kx, ky, 
                    berry_s[0].berry,
                    berry_s[1].berry,
                    berry_s[2].berry,
                    berry_s[3].berry,
                    berry_s[4].berry,
                    berry_s[5].berry,
                    berry_s[0].eigen,
                    berry_s[1].eigen,
                    berry_s[2].eigen,
                    berry_s[3].eigen,
                    berry_s[4].eigen,
                    berry_s[5].eigen,
                    berry_s[0].spin.to_u8(),
                    berry_s[1].spin.to_u8(),
                    berry_s[2].spin.to_u8(),
                    berry_s[3].spin.to_u8(),
                    berry_s[4].spin.to_u8(),
                    berry_s[5].spin.to_u8(),
                ).unwrap();
            }
        

        }
    }

    if export_dat{
        fs::write(
            format!(
                "./output/band_infos/dats/{}_{}_lambda_{}_j{}_mu{}.dat"
                ,comment
                ,graph_mesh
                ,lambda.to_string().replace('.', "p")
                ,jj.to_string().replace('.', "p")
                ,mu.to_string().replace('.', "p")
            ), 
            file_str
        ).unwrap();
    }

}


pub fn point_in_triangle_simple(
    x: Vector2<f64>,
    a: Vector2<f64>,
    b: Vector2<f64>,
    c: Vector2<f64>,
) -> bool {
    const EPSILON: f64 = 0.;

    let cross_ab_ax = (b.x - a.x) * (x.y - a.y) - (b.y - a.y) * (x.x - a.x);
    let cross_bc_bx = (c.x - b.x) * (x.y - b.y) - (c.y - b.y) * (x.x - b.x);
    let cross_ca_cx = (a.x - c.x) * (x.y - c.y) - (a.y - c.y) * (x.x - c.x);

    let all_non_negative = cross_ab_ax >= -EPSILON 
                        && cross_bc_bx >= -EPSILON 
                        && cross_ca_cx >= -EPSILON;

    let all_non_positive = cross_ab_ax <= EPSILON 
                        && cross_bc_bx <= EPSILON 
                        && cross_ca_cx <= EPSILON;

    all_non_negative || all_non_positive
}