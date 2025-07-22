use crate::consts::*;
use crate::diag::{diag_66,diag_22,SEud6,SEud2};
use crate::util::point_in_triangle_simple;
use crate::dv2::DV2;
use crate::model::System;

use nalgebra::{Complex, Vector2, Vector6};
use std::fs;
use std::fmt::Write;


#[derive(Debug, Clone, Copy)]
struct Binfo{
    berry : f64,
    eigen : f64,
    spin : Spin,
    bcd : Vector2<f64>
}

impl Binfo{
    fn ini() -> Self{
        Binfo { berry: 0.0, eigen: 0.0, spin : Spin::Undefined, bcd : Vector2::zeros() }
    }
    fn new(berry : f64, eigen : f64, spin : Spin, bcd : Vector2<f64>) -> Self{
        Binfo { berry, eigen , spin, bcd}
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

pub fn calculate_band_info_all_band(export_dat : bool, graph_mesh : usize, system : &System, comment : &str){

    let mesh_kx = graph_mesh;
    let mesh_ky  = graph_mesh;

    match system.size(){
        6 => {
            let dv2_k1 = DV2::from_car(KPPKS) - DV2::from_car(-KINKS);
            let dv2_k2 = DV2::from_car(KP_KS) - DV2::from_car(-KINKS);

            let mut seud_mat = vec![vec![None ; mesh_ky + 2] ; mesh_kx + 2];

            for i in 0..(mesh_kx+2){
                for j in 0..(mesh_ky+2){
                    let if64 = i as f64 / mesh_kx as f64;
                    let jf64 = j as f64 / mesh_ky as f64;

                    let kk_dv2 = dv2_k1 * if64 + dv2_k2 * jf64 + DV2::from_car(-KINKS);

                    let seud = diag_66(kk_dv2.to_car(),system).sort();

                    seud_mat[i][j] = Some(seud);
                }
            }

            let mut file_str = format!("kx,ky,b1,b2,b3,b4,b5,b6,e1,e2,e3,e4,e5,s1,s2,s3,s4,s5,s6,bcd1x,bcd1y,bcd2x,bcd2y,bcd3x,bcd3y,bcd4x,bcd4y,bcd5x,bcd5y,bcd6x,bcd6y \n{} {}",graph_mesh,system.debug());

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

                    let lattice_1_len = dv2_k1.to_car().norm() / mesh_kx as f64;
                    let lattice_2_len = dv2_k2.to_car().norm() / mesh_ky as f64;

                    let cell_area = lattice_1_len * lattice_2_len * SQRT_3 * 0.5;
                    

                    let kx = kk_dv2.to_car().x;
                    let ky = kk_dv2.to_car().y;

                    let berrys = seud_mat_to_bc6(i, j, &seud_mat);
                    let berrys_pi = seud_mat_to_bc6(i + 1, j, &seud_mat);
                    let berrys_pi2 = seud_mat_to_bc6(i + 1, (j + mesh_ky - 1) % mesh_ky, &seud_mat);
                    let berrys_pj = seud_mat_to_bc6(i , j + 1, &seud_mat);

                    let seud = SEud6::get_cont(&seud_mat[i][j]);

                    let se_u = seud.u;
                    let se_d = seud.d;

                    let mut binfos = [Binfo::ini();12];

                    for eigen_i in 0..se_u.eigenvalues.len(){
                        {
                            let berry = berrys.up_berrys[eigen_i] / cell_area;

                            let bcd = {
                                
                                let x = (berrys_pj.up_berrys[eigen_i] - berrys.up_berrys[eigen_i]) / lattice_2_len;

                                let ave = (berrys_pi.up_berrys[eigen_i] + berrys_pi2.up_berrys[eigen_i]) * 0.5;
                                let y = ((berrys.up_berrys[eigen_i] - ave) / lattice_1_len) * 2. / SQRT_3;

                                Vector2::new(x,y)
                            };

                            binfos[eigen_i] = Binfo::new(berry,se_u.eigenvalues[eigen_i],Spin::U,bcd);
                        }
                    }

                    for eigen_i in 0..se_d.eigenvalues.len(){
                        {
                            let berry = berrys.do_berrys[eigen_i] / cell_area;

                            let bcd = {
                                
                                let x = (berrys_pj.do_berrys[eigen_i] - berrys.do_berrys[eigen_i]) / lattice_2_len;

                                let ave = (berrys_pi.do_berrys[eigen_i] + berrys_pi2.do_berrys[eigen_i]) * 0.5;
                                let y = ((berrys.do_berrys[eigen_i] - ave) / lattice_1_len) * 2. / SQRT_3;

                                Vector2::new(x,y)
                            };

                            binfos[eigen_i + 6] = Binfo::new(berry,se_d.eigenvalues[eigen_i],Spin::D,bcd);
                        }
                    }

                    binfos.sort_by(|a, b| a.eigen.partial_cmp(&b.eigen).unwrap());


                    if export_dat{
                        write!(file_str, "\n{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}", 
                            kx, ky, 
                            binfos[0].berry,
                            binfos[1].berry,
                            binfos[2].berry,
                            binfos[3].berry,
                            binfos[4].berry,
                            binfos[5].berry,
                            binfos[0].eigen,
                            binfos[1].eigen,
                            binfos[2].eigen,
                            binfos[3].eigen,
                            binfos[4].eigen,
                            binfos[5].eigen,
                            binfos[0].spin.to_u8(),
                            binfos[1].spin.to_u8(),
                            binfos[2].spin.to_u8(),
                            binfos[3].spin.to_u8(),
                            binfos[4].spin.to_u8(),
                            binfos[5].spin.to_u8(),
                            binfos[0].bcd.x,
                            binfos[0].bcd.y,
                            binfos[1].bcd.x,
                            binfos[1].bcd.y,
                            binfos[2].bcd.x,
                            binfos[2].bcd.y,
                            binfos[3].bcd.x,
                            binfos[3].bcd.y,
                            binfos[4].bcd.x,
                            binfos[4].bcd.y,
                            binfos[5].bcd.x,
                            binfos[5].bcd.y,
                        ).unwrap();
                    }
                

                }
            }

            if export_dat{
                fs::write(
                    format!(
                        "./output/band_infos/dats/{}_{}_{}.csv"
                        ,comment
                        ,graph_mesh
                        ,system.debug()
                    ), 
                    file_str
                ).unwrap();
            }
        }
        2 => {
            let dv2_k1 = DV2::from_car2(KPPKS2) - DV2::from_car2(-KINKS2);
            let dv2_k2 = DV2::from_car2(KP_KS2) - DV2::from_car2(-KINKS2);

            let mut seud_mat = vec![vec![None ; mesh_ky + 2] ; mesh_kx + 2];

            for i in 0..(mesh_kx+2){
                for j in 0..(mesh_ky+2){
                    let if64 = i as f64 / mesh_kx as f64;
                    let jf64 = j as f64 / mesh_ky as f64;

                    let kk_dv2 = dv2_k1 * if64 + dv2_k2 * jf64 + DV2::from_car2(-KINKS2);

                    let seud = diag_22(kk_dv2.to_car2(),system).sort();

                    seud_mat[i][j] = Some(seud);
                }
            }

            let mut file_str = format!("kx,ky,b1,b2,e1,e2,s1,s2,bcd1x,bcd1y,bcd2x,bcd2y\n{} {}",graph_mesh,system.debug());

            for i in 0..mesh_kx{
                for j in 0..mesh_ky{
                    let if64 = i as f64 / mesh_kx as f64;
                    let jf64 = j as f64 / mesh_ky as f64;

                    let mut kk_dv2 = dv2_k1 * if64 + dv2_k2 * jf64 + DV2::from_car2(-KINKS2);

                    if point_in_triangle_simple(kk_dv2.to_car2(), KINKS2, 2. * KINKS2, KPPKS2){
                        kk_dv2 = kk_dv2 - dv2_k1;
                    }
                    else if point_in_triangle_simple(kk_dv2.to_car2(), KINKS2, 2. * KINKS2, KP_KS2){
                        kk_dv2 = kk_dv2 - dv2_k2;
                    }

                    let lattice_1_len = dv2_k1.to_car2().norm() / mesh_kx as f64;
                    let lattice_2_len = dv2_k2.to_car2().norm() / mesh_ky as f64;

                    let cell_area = lattice_1_len * lattice_2_len * SQRT_3 * 0.5;
                    

                    let kx = kk_dv2.to_car2().x;
                    let ky = kk_dv2.to_car2().y;

                    let berrys = seud_mat_to_bc2(i, j, &seud_mat);
                    let berrys_pi = seud_mat_to_bc2(i + 1, j, &seud_mat);
                    let berrys_pimj = seud_mat_to_bc2(i + 1, (j - 1 + mesh_ky) % mesh_ky, &seud_mat);
                    let berrys_pj = seud_mat_to_bc2(i , j + 1, &seud_mat);

                    let seud = SEud2::get_cont(&seud_mat[i][j]);

                    let se_u = seud.u;
                    let se_d = seud.d;

                    let mut binfos = [Binfo::ini();4];

                    for eigen_i in 0..se_u.eigenvalues.len(){
                        {
                            let berry = berrys.up_berrys[eigen_i] / cell_area;

                            let bcd = {
                                let y = (berrys.up_berrys[eigen_i] - berrys_pimj.up_berrys[eigen_i]) / lattice_1_len;

                                let ave = (berrys_pi.up_berrys[eigen_i] + berrys_pj.up_berrys[eigen_i]) * 0.5;
                                let x = ((ave - berrys.up_berrys[eigen_i]) / lattice_2_len) * 2. / SQRT_3;

                                Vector2::new(x,y)
                            };

                            binfos[eigen_i] = Binfo::new(berry,se_u.eigenvalues[eigen_i],Spin::U,bcd);
                        }
                    }

                    for eigen_i in 0..se_d.eigenvalues.len(){
                        {
                            let berry = berrys.do_berrys[eigen_i] / cell_area;

                            let bcd = {
                                let y = (berrys.do_berrys[eigen_i] - berrys_pimj.do_berrys[eigen_i]) / lattice_1_len;

                                let ave = (berrys_pi.do_berrys[eigen_i] + berrys_pj.do_berrys[eigen_i]) * 0.5;
                                let x = ((ave - berrys.do_berrys[eigen_i]) / lattice_2_len) * 2. / SQRT_3;

                                Vector2::new(x,y)
                            };

                            binfos[eigen_i + 2] = Binfo::new(berry,se_d.eigenvalues[eigen_i],Spin::D,bcd);
                        }
                    }

                    binfos.sort_by(|a, b| a.eigen.partial_cmp(&b.eigen).unwrap());


                    if export_dat{
                        write!(file_str, "\n{},{},{},{},{},{},{},{},{},{},{},{}", 
                            kx, ky, 
                            binfos[0].berry,
                            binfos[1].berry,
                            binfos[0].eigen,
                            binfos[1].eigen,
                            binfos[0].spin.to_u8(),
                            binfos[1].spin.to_u8(),
                            binfos[0].bcd.x,
                            binfos[0].bcd.y,
                            binfos[1].bcd.x,
                            binfos[1].bcd.y,
                        ).unwrap();
                    }
                

                }
            }

            if export_dat{
                fs::write(
                    format!(
                        "./output/band_infos/dats/{}_{}_{}.csv"
                        ,comment
                        ,graph_mesh
                        ,system.debug()
                    ), 
                    file_str
                ).unwrap();
            }
        }
        _ => {}
    }
    
    

}

struct Berrys6{
    up_berrys : [f64;6],
    do_berrys : [f64;6]
}

fn seud_mat_to_bc6(i : usize,j : usize ,seud_mat : &Vec<Vec<Option<SEud6>>>) -> Berrys6{
    let seud = SEud6::get_cont(&seud_mat[i][j]);

    let se_u = seud.u;
    let se_d = seud.d;

    let mut berrys = Berrys6{up_berrys : [0.0;6],do_berrys : [0.0;6]};

    for eigen_i in 0..se_u.eigenvalues.len(){
        {
            let u00 : Vector6<Complex<f64>> = se_u.eigenvectors.column(eigen_i).into();
            let u10 : Vector6<Complex<f64>> = SEud6::get_cont(&seud_mat[i + 1][j]).u.eigenvectors.column(eigen_i).into();
            let u11 : Vector6<Complex<f64>> = SEud6::get_cont(&seud_mat[i + 1][j + 1]).u.eigenvectors.column(eigen_i).into();
            let u01 : Vector6<Complex<f64>> = SEud6::get_cont(&seud_mat[i][j + 1]).u.eigenvectors.column(eigen_i).into();

            let berry = {
                let u1 = u00.dotc(&u10);
                let u2 = u10.dotc(&u11);
                let u3 = u11.dotc(&u01);
                let u4 = u01.dotc(&u00);

                (u1 * u2 * u3 * u4).arg()
            };

            berrys.up_berrys[eigen_i] = berry;
        }
    }

    for eigen_i in 0..se_d.eigenvalues.len(){
        {
            let u00 : Vector6<Complex<f64>> = se_d.eigenvectors.column(eigen_i).into();
            let u10 : Vector6<Complex<f64>> = SEud6::get_cont(&seud_mat[i + 1][j]).d.eigenvectors.column(eigen_i).into();
            let u11 : Vector6<Complex<f64>> = SEud6::get_cont(&seud_mat[i + 1][j + 1]).d.eigenvectors.column(eigen_i).into();
            let u01 : Vector6<Complex<f64>> = SEud6::get_cont(&seud_mat[i][j + 1]).d.eigenvectors.column(eigen_i).into();

            let berry = {
                let u1 = u00.dotc(&u10);
                let u2 = u10.dotc(&u11);
                let u3 = u11.dotc(&u01);
                let u4 = u01.dotc(&u00);

                (u1 * u2 * u3 * u4).arg()
            };

            berrys.do_berrys[eigen_i] = berry;
        }
    }

    berrys
}

struct Berrys2{
    up_berrys : [f64;6],
    do_berrys : [f64;6]
}

fn seud_mat_to_bc2(i : usize,j : usize ,seud_mat : &Vec<Vec<Option<SEud2>>>) -> Berrys2{
    let seud = SEud2::get_cont(&seud_mat[i][j]);

    let se_u = seud.u;
    let se_d = seud.d;

    let mut berrys = Berrys2{up_berrys : [0.0;6],do_berrys : [0.0;6]};

    for eigen_i in 0..se_u.eigenvalues.len(){
        {
            let u00 : Vector2<Complex<f64>> = se_u.eigenvectors.column(eigen_i).into();
            let u10 : Vector2<Complex<f64>> = SEud2::get_cont(&seud_mat[i + 1][j]).u.eigenvectors.column(eigen_i).into();
            let u11 : Vector2<Complex<f64>> = SEud2::get_cont(&seud_mat[i + 1][j + 1]).u.eigenvectors.column(eigen_i).into();
            let u01 : Vector2<Complex<f64>> = SEud2::get_cont(&seud_mat[i][j + 1]).u.eigenvectors.column(eigen_i).into();

            let berry = {
                let u1 = u00.dotc(&u10);
                let u2 = u10.dotc(&u11);
                let u3 = u11.dotc(&u01);
                let u4 = u01.dotc(&u00);

                (u1 * u2 * u3 * u4).arg()
            };

            berrys.up_berrys[eigen_i] = berry;
        }
    }

    for eigen_i in 0..se_d.eigenvalues.len(){
        {
            let u00 : Vector2<Complex<f64>> = se_d.eigenvectors.column(eigen_i).into();
            let u10 : Vector2<Complex<f64>> = SEud2::get_cont(&seud_mat[i + 1][j]).d.eigenvectors.column(eigen_i).into();
            let u11 : Vector2<Complex<f64>> = SEud2::get_cont(&seud_mat[i + 1][j + 1]).d.eigenvectors.column(eigen_i).into();
            let u01 : Vector2<Complex<f64>> = SEud2::get_cont(&seud_mat[i][j + 1]).d.eigenvectors.column(eigen_i).into();

            let berry = {
                let u1 = u00.dotc(&u10);
                let u2 = u10.dotc(&u11);
                let u3 = u11.dotc(&u01);
                let u4 = u01.dotc(&u00);

                (u1 * u2 * u3 * u4).arg()
            };

            berrys.do_berrys[eigen_i] = berry;
        }
    }

    berrys
}