use crate::consts::*;
use crate::diag::{diag,SEudEnum};
use crate::util::point_in_triangle_simple;
use crate::dv2::DV2;
use crate::model::System;

use nalgebra::{Complex, Vector2, Vector6,};
use std::fs;
use std::fmt::Write;

const DELTA : f64 = 0.001;

struct Binfos{
    u : Vec<Binfo>,
    d : Vec<Binfo>
}
impl Binfos {
    pub fn merged(&self) -> Vec<Binfo> {
        [self.u.clone(), self.d.clone()].concat()
    }
}


#[derive(Debug, Clone, Copy)]
struct Binfo{
    berry : f64,
    eigen : f64,
    spin : Spin,
    bcd : Option<Vector2<f64>>
}

impl Binfo{
    fn ini() -> Self{
        Binfo { berry: 0.0, eigen: 0.0, spin : Spin::Undefined, bcd : None }
    }
    fn new(berry : f64, eigen : f64, spin : Spin, bcd : Option<Vector2<f64>>) -> Self{
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

pub fn calculate_band_info_all_band(graph_mesh : usize, system : &System, comment : &str){

    let mut file_str = format!("{}\n{} {}\n"
        ,{
            match system.size(){
                2 => {"kx,ky,b1,b2,e1,e2,s1,s2,bcd1x,bcd1y,bcd2x,bcd2y"}
                6 => {"x,ky,b1,b2,b3,b4,b5,b6,e1,e2,e3,e4,e5,s1,s2,s3,s4,s5,s6,bcd1x,bcd1y,bcd2x,bcd2y,bcd3x,bcd3y,bcd4x,bcd4y,bcd5x,bcd5y,bcd6x,bcd6y"}
                _ => {panic!("invalid sizes")}
            }
        }
        ,graph_mesh,system.debug()
    );

    let mesh_kx = graph_mesh;
    let mesh_ky  = graph_mesh;

    let mut seud_mat = vec![vec![None ; mesh_ky + 1] ; mesh_kx + 1];
    let mut seud_mat_px = vec![vec![None ; mesh_ky + 1] ; mesh_kx + 1];
    let mut seud_mat_py = vec![vec![None ; mesh_ky + 1] ; mesh_kx + 1];

    for i in 0..(mesh_kx+1){
        for j in 0..(mesh_ky+1){
            let kk = i_j_to_kk(i, j, mesh_kx, mesh_ky, false, system.size()).0;

            let seud = diag(kk,system).sort();
            let seud_px = diag(kk + Vector2::new(DELTA,0.0),system).sort();
            let seud_py = diag(kk + Vector2::new(0.0,DELTA),system).sort();

            seud_mat[i][j] = Some(seud);
            seud_mat_px[i][j] = Some(seud_px);
            seud_mat_py[i][j] = Some(seud_py);
        }
    }

    for i in 0..(mesh_kx){
        for j in 0..(mesh_ky){
            let kk_and_cell_area = i_j_to_kk(i, j, mesh_kx, mesh_ky, true, system.size());
            let kk = kk_and_cell_area.0;
            
            let band_infos = cal_band_infos(i,j,&seud_mat,&seud_mat_px,&seud_mat_py,system.size(),kk_and_cell_area.1);

            let berrys : Vec<f64> = band_infos.iter().map(|b| b.berry).collect();
            let eigens : Vec<f64> = band_infos.iter().map(|b| b.eigen).collect();
            let spins  : Vec<u8> = band_infos.iter().map(|b| b.spin.to_u8()).collect();
            let bcds   : Vec<Vector2<f64>> = band_infos.iter().map(|b| b.bcd.unwrap()).collect();

            write!(file_str,"{},{}",kk.x,kk.y).unwrap();

            for k in 0..system.size(){
                write!(file_str,",{}",
                    berrys[k],
                ).unwrap();
            }
            for k in 0..system.size(){
                write!(file_str,",{}",
                    eigens[k],
                ).unwrap();
            }
            for k in 0..system.size(){
                write!(file_str,",{}",
                    spins[k],
                ).unwrap();
            }
            for k in 0..system.size(){
                write!(file_str,",{},{}",
                    bcds[k].x,bcds[k].y
                ).unwrap();
            }

            write!(file_str,"\n").unwrap();
        }
    }

    fs::write(
        format!(
            "./output/band_infos/dats/{}_{}_{}.csv"
            ,comment
            ,graph_mesh
            ,system.debug()
        )
        , file_str
    ).unwrap();

}


fn cal_band_infos(
    i : usize, 
    j : usize, 
    seud_mat : &Vec<Vec<Option<SEudEnum>>>,
    seud_mat_px : &Vec<Vec<Option<SEudEnum>>>,
    seud_mat_py : &Vec<Vec<Option<SEudEnum>>>,
    size : usize,
    cell_area : f64,
) -> Vec<Binfo>{
    let mut berrys = cal_bes(i, j, seud_mat,size,cell_area);

    let berrys_px = cal_bes(i, j, seud_mat_px,size,cell_area);
    let berrys_py = cal_bes(i, j, seud_mat_py,size,cell_area);

    for ei in 0..size{
        let berry = berrys.u[ei].berry;
        let berry_px = berrys_px.u[ei].berry;
        let berry_py = berrys_py.u[ei].berry;

        berrys.u[ei].bcd = Some({
            let changex = (berry_px - berry) / DELTA;
            let changey = (berry_py - berry) / DELTA;

            Vector2::new(changex,changey)
        })
    }
    for ei in 0..size{
        let berry = berrys.d[ei].berry;
        let berry_px = berrys_px.d[ei].berry;
        let berry_py = berrys_py.d[ei].berry;

        berrys.d[ei].bcd = Some({
            let changex = (berry_px - berry) / DELTA;
            let changey = (berry_py - berry) / DELTA;

            Vector2::new(changex,changey)
        })
    }

    let mut vec_berrys = berrys.merged();

    vec_berrys.sort_by(|a, b| a.eigen.partial_cmp(&b.eigen).unwrap());

    vec_berrys
}

fn cal_bes(
    i : usize, 
    j : usize, 
    seud_mat : &Vec<Vec<Option<SEudEnum>>>,
    size : usize,
    cell_area : f64,
) -> Binfos {

    match size{
        2 => {
            let seud = SEudEnum::get_cont(&seud_mat[i][j]).is_2();
            let seud_pi = SEudEnum::get_cont(&seud_mat[i + 1][j]).is_2();
            let seud_pij = SEudEnum::get_cont(&seud_mat[i + 1][j + 1]).is_2();
            let seud_pj = SEudEnum::get_cont(&seud_mat[i][j + 1]).is_2();

            let mut berrys_up = vec![Binfo::ini();2];
            let mut berrys_do = vec![Binfo::ini();2];

            for ei in 0..2{
                let u00 : Vector2<Complex<f64>> = seud.u.eigenvectors.column(ei).into();
                let u10 : Vector2<Complex<f64>> = seud_pi.u.eigenvectors.column(ei).into();
                let u11 : Vector2<Complex<f64>> = seud_pij.u.eigenvectors.column(ei).into();
                let u01 : Vector2<Complex<f64>> = seud_pj.u.eigenvectors.column(ei).into(); 

                let berry = {
                    let u1 = u00.dotc(&u10);
                    let u2 = u10.dotc(&u11);
                    let u3 = u11.dotc(&u01);
                    let u4 = u01.dotc(&u00);

                    (u1 * u2 * u3 * u4).arg()
                } / cell_area;

                berrys_up[ei] = Binfo::new(berry, seud.u.eigenvalues[ei],Spin::U,None)

            }
            for ei in 0..2{
                let u00 : Vector2<Complex<f64>> = seud.d.eigenvectors.column(ei).into();
                let u10 : Vector2<Complex<f64>> = seud_pi.d.eigenvectors.column(ei).into();
                let u11 : Vector2<Complex<f64>> = seud_pij.d.eigenvectors.column(ei).into();
                let u01 : Vector2<Complex<f64>> = seud_pj.d.eigenvectors.column(ei).into(); 

                let berry = {
                    let u1 = u00.dotc(&u10);
                    let u2 = u10.dotc(&u11);
                    let u3 = u11.dotc(&u01);
                    let u4 = u01.dotc(&u00);

                    (u1 * u2 * u3 * u4).arg()
                } / cell_area;

                berrys_do[ei] = Binfo::new(berry, seud.d.eigenvalues[ei],Spin::D,None)

            }

            Binfos{u : berrys_up, d : berrys_do}
        }
        6 => {
            let seud = SEudEnum::get_cont(&seud_mat[i][j]).is_6();
            let seud_pi = SEudEnum::get_cont(&seud_mat[i + 1][j]).is_6();
            let seud_pij = SEudEnum::get_cont(&seud_mat[i + 1][j + 1]).is_6();
            let seud_pj = SEudEnum::get_cont(&seud_mat[i][j + 1]).is_6();

            let mut berrys_up = vec![Binfo::ini();6];
            let mut berrys_do = vec![Binfo::ini();6];

            for ei in 0..6{
                let u00 : Vector6<Complex<f64>> = seud.u.eigenvectors.column(ei).into();
                let u10 : Vector6<Complex<f64>> = seud_pi.u.eigenvectors.column(ei).into();
                let u11 : Vector6<Complex<f64>> = seud_pij.u.eigenvectors.column(ei).into();
                let u01 : Vector6<Complex<f64>> = seud_pj.u.eigenvectors.column(ei).into(); 

                let berry = {
                    let u1 = u00.dotc(&u10);
                    let u2 = u10.dotc(&u11);
                    let u3 = u11.dotc(&u01);
                    let u4 = u01.dotc(&u00);

                    (u1 * u2 * u3 * u4).arg()
                } / cell_area;

                berrys_up[ei] = Binfo::new(berry, seud.u.eigenvalues[ei],Spin::U,None)

            }
            for ei in 0..6{
                let u00 : Vector6<Complex<f64>> = seud.d.eigenvectors.column(ei).into();
                let u10 : Vector6<Complex<f64>> = seud_pi.d.eigenvectors.column(ei).into();
                let u11 : Vector6<Complex<f64>> = seud_pij.d.eigenvectors.column(ei).into();
                let u01 : Vector6<Complex<f64>> = seud_pj.d.eigenvectors.column(ei).into(); 

                let berry = {
                    let u1 = u00.dotc(&u10);
                    let u2 = u10.dotc(&u11);
                    let u3 = u11.dotc(&u01);
                    let u4 = u01.dotc(&u00);

                    (u1 * u2 * u3 * u4).arg()
                } / cell_area;

                berrys_do[ei] = Binfo::new(berry, seud.d.eigenvalues[ei],Spin::D,None)

            }

            Binfos{u : berrys_up, d : berrys_do}
        }
        _ => panic!("invalid size")
    }
}

fn i_j_to_kk(i : usize, j : usize, mesh_kx : usize, mesh_ky : usize, hex : bool, size : usize) -> (Vector2<f64>,f64){
    match size{
        6 => {
            let dv2_k1 = DV2::from_car(KPPKS) - DV2::from_car(-KINKS);
            let dv2_k2 = DV2::from_car(KP_KS) - DV2::from_car(-KINKS);

            let if64 = i as f64 / mesh_kx as f64;
            let jf64 = j as f64 / mesh_ky as f64;

            let mut kk_dv2 = dv2_k1 * if64 + dv2_k2 * jf64 + DV2::from_car(-KINKS);

            if hex{
                if point_in_triangle_simple(kk_dv2.to_car(), KINKS, 2. * KINKS, KPPKS){
                    kk_dv2 = kk_dv2 - dv2_k1;
                }
                else if point_in_triangle_simple(kk_dv2.to_car(), KINKS, 2. * KINKS, KP_KS){
                    kk_dv2 = kk_dv2 - dv2_k2;
                }
            }

            let lattice_1_len = dv2_k1.to_car().norm() / mesh_kx as f64;
            let lattice_2_len = dv2_k2.to_car().norm() / mesh_ky as f64;

            let cell_area = lattice_1_len * lattice_2_len * SQRT_3 * 0.5;
            

            let kx = kk_dv2.to_car().x;
            let ky = kk_dv2.to_car().y;

            (Vector2::new(kx,ky),cell_area)
        }
        2 => {
            let dv2_k1 = DV2::from_car2(KPPKS2) - DV2::from_car2(-KINKS2);
            let dv2_k2 = DV2::from_car2(KP_KS2) - DV2::from_car2(-KINKS2);

            let if64 = i as f64 / mesh_kx as f64;
            let jf64 = j as f64 / mesh_ky as f64;

            let mut kk_dv2 = dv2_k1 * if64 + dv2_k2 * jf64 + DV2::from_car2(-KINKS2);

            if hex{
                if point_in_triangle_simple(kk_dv2.to_car2(), KINKS2, 2. * KINKS2, KPPKS2){
                    kk_dv2 = kk_dv2 - dv2_k1;
                }
                else if point_in_triangle_simple(kk_dv2.to_car2(), KINKS2, 2. * KINKS2, KP_KS2){
                    kk_dv2 = kk_dv2 - dv2_k2;
                }
            }

            let lattice_1_len = dv2_k1.to_car2().norm() / mesh_kx as f64;
            let lattice_2_len = dv2_k2.to_car2().norm() / mesh_ky as f64;

            let cell_area = lattice_1_len * lattice_2_len * SQRT_3 * 0.5;
            

            let kx = kk_dv2.to_car2().x;
            let ky = kk_dv2.to_car2().y;

            (Vector2::new(kx,ky),cell_area)
        }
        _ => {panic!("invalid size")}
    }
}