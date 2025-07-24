use crate::consts::*;
use crate::diag::{diag,SEudEnum};
use crate::util::{i_j_to_kk,cal_cell_area,cal_cell_area_grid,i_j_to_kk_grid,GridInfo};
use crate::model::System;
use crate::binfo::BinfosMergedOnkk;

use crate::cal_bc_bcd::{bc_sum_up,bcd_sum_up};

use nalgebra::{Vector2,};
use std::fs;

pub fn cal_bc_grid_test(fermi_energy : f64, system : &System) -> f64{

    let mut bc = 0.;

    let n = 100;

    for grid_x in 0..n{
        for grid_y in 0..n{
            let grid_info = GridInfo::new_ijn(grid_x,grid_y,n,n);

            let (binfos_merged_onkks, max_bc ) = calculate_band_info_grid(3, system, grid_info);

            if max_bc.abs() > 100.0{
                let (binfos_merged_onkks_new, _ ) = calculate_band_info_grid(1000, system, grid_info);

                bc += bc_sum_up(fermi_energy, &binfos_merged_onkks_new);
            }
            else if max_bc.abs() > 10.0{
                let (binfos_merged_onkks_new, _ ) = calculate_band_info_grid(100, system, grid_info);

                bc += bc_sum_up(fermi_energy, &binfos_merged_onkks_new);
            }
            else if max_bc.abs() > 1.0{
                let (binfos_merged_onkks_new, _ ) = calculate_band_info_grid(10, system, grid_info);

                bc += bc_sum_up(fermi_energy, &binfos_merged_onkks_new);
            }
            else{
                bc += bc_sum_up(fermi_energy, &binfos_merged_onkks);
            }
        }
    }

    bc
}

pub fn cal_bcd_grid_test(fermi_energy : f64, system : &System) -> Vector2<f64>{

    let mut bcd = Vector2::zeros();

    let n = 500;

    for grid_x in 0..n{
        for grid_y in 0..n{
            let grid_info = GridInfo::new_ijn(grid_x,grid_y,n,n);

            let (binfos_merged_onkks, max_bcd ) = calculate_band_info_grid(1, system, grid_info);

            if max_bcd.abs() > 10.{

                let graph_mesh = std::cmp::min((max_bcd.abs()/10.) as usize,1000);
                let (binfos_merged_onkks_new, _ ) = calculate_band_info_grid(graph_mesh, system, grid_info);

                bcd += bcd_sum_up(fermi_energy, &binfos_merged_onkks_new);
            }
            else{
                bcd += bcd_sum_up(fermi_energy, &binfos_merged_onkks);
            }
        }
    }

    bcd
}

pub fn cal_bcd_grid(fermi_energy : f64, system : &System,grid_mesh : usize, max_sub_mesh : usize, mesh_scale : f64) -> Vector2<f64>{

    let mut bcd = Vector2::zeros();

    let n = grid_mesh;

    for grid_x in 0..n{
        for grid_y in 0..n{
            let grid_info = GridInfo::new_ijn(grid_x,grid_y,n,n);

            let (binfos_merged_onkks, max_bcd ) = calculate_band_info_grid(1, system, grid_info);

            if max_bcd.abs() > 10.{

                let graph_mesh = std::cmp::min((max_bcd.abs() * mesh_scale) as usize,max_sub_mesh);
                let (binfos_merged_onkks_new, _ ) = calculate_band_info_grid(graph_mesh, system, grid_info);

                bcd += bcd_sum_up(fermi_energy, &binfos_merged_onkks_new);
            }
            else{
                bcd += bcd_sum_up(fermi_energy, &binfos_merged_onkks);
            }
        }
    }

    bcd
}

pub fn calculate_band_info_grid(graph_mesh : usize, system : &System, grid_info : GridInfo) -> (Vec<BinfosMergedOnkk>,f64) {
    let mut out = Vec::new();
    let mut max_bc = 0.0f64;
    let mut max_bcd = 0.0f64;

    let mesh_kx = graph_mesh;
    let mesh_ky  = graph_mesh;

    let seud_mat = make_seud_mat_grid(Vector2::zeros(), mesh_kx, mesh_ky, system,grid_info);
    let seud_mat_px = make_seud_mat_grid(Vector2::new(DELTA,0.0), mesh_kx, mesh_ky, system,grid_info);
    let seud_mat_py = make_seud_mat_grid(Vector2::new(0.0,DELTA), mesh_kx, mesh_ky, system,grid_info);
    let seud_mat_mx = make_seud_mat_grid(Vector2::new(-DELTA,0.0), mesh_kx, mesh_ky, system,grid_info);
    let seud_mat_my = make_seud_mat_grid(Vector2::new(0.0,-DELTA), mesh_kx, mesh_ky, system,grid_info);

    for i in 0..(mesh_kx){
        for j in 0..(mesh_ky){
            let kk = i_j_to_kk_grid(i, j, mesh_kx, mesh_ky, true, system.size(),grid_info);
            
            let binfos_merged_onkk = BinfosMergedOnkk::cal_cd(i,j,&seud_mat,&seud_mat_px,&seud_mat_py,&seud_mat_mx,&seud_mat_my,system.size(),kk);

            let cell_area = cal_cell_area_grid(mesh_kx, mesh_ky, system.size(), grid_info);

            let max_bc_in_this_kk = binfos_merged_onkk.max_bc(cell_area);
            let max_bcd_in_this_kk = binfos_merged_onkk.max_bcd(cell_area);

            if max_bc.abs() < max_bc_in_this_kk.abs() {
                max_bc = max_bc_in_this_kk;
            }
            if max_bcd.abs() < max_bcd_in_this_kk.abs() {
                max_bcd = max_bcd_in_this_kk;
            }

            out.push(binfos_merged_onkk);
        }
    }

    println!("{}",max_bcd);

    (out,max_bcd)
}

//BZを等間隔に切ってBCやBCDを計算
pub fn calculate_band_info_all_band(write : bool ,graph_mesh : usize, system : &System, comment : &str) -> Vec<BinfosMergedOnkk>{

    let mut out = Vec::new();

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

    let seud_mat = make_seud_mat(Vector2::zeros(), mesh_kx, mesh_ky, system);
    let seud_mat_px = make_seud_mat(Vector2::new(DELTA,0.0), mesh_kx, mesh_ky, system);
    let seud_mat_py = make_seud_mat(Vector2::new(0.0,DELTA), mesh_kx, mesh_ky, system);
    let seud_mat_mx = make_seud_mat(Vector2::new(-DELTA,0.0), mesh_kx, mesh_ky, system);
    let seud_mat_my = make_seud_mat(Vector2::new(0.0,-DELTA), mesh_kx, mesh_ky, system);

    let cell_area = cal_cell_area(mesh_kx, mesh_ky, system.size());

    for i in 0..(mesh_kx){
        for j in 0..(mesh_ky){
            let kk = i_j_to_kk(i, j, mesh_kx, mesh_ky, true, system.size());
            
            let binfos_merged_onkk = BinfosMergedOnkk::cal_cd(i,j,&seud_mat,&seud_mat_px,&seud_mat_py,&seud_mat_mx,&seud_mat_my,system.size(),kk);

            if write{
                binfos_merged_onkk.write(&mut file_str, system.size(), cell_area);
            }

            out.push(binfos_merged_onkk);
        }
    }

    if write{
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

    out

}

fn make_seud_mat_grid(deltakk : Vector2<f64>, mesh_kx : usize, mesh_ky : usize,system : &System, grid_info : GridInfo) -> Vec<Vec<Option<SEudEnum>>>{
    let mut seud_mat = vec![vec![None ; mesh_ky + 1] ; mesh_kx + 1];

    for i in 0..(mesh_kx + 1){
        for j in 0..(mesh_ky + 1){
            let kk = i_j_to_kk_grid(i, j, mesh_kx, mesh_ky, false, system.size(),grid_info);

            let seud = diag(kk + deltakk,system).sort();

            seud_mat[i][j] = Some(seud);
        }
    }

    seud_mat
}

fn make_seud_mat(deltakk : Vector2<f64>, mesh_kx : usize, mesh_ky : usize,system : &System) -> Vec<Vec<Option<SEudEnum>>>{
    let mut seud_mat = vec![vec![None ; mesh_ky + 1] ; mesh_kx + 1];

    for i in 0..(mesh_kx + 1){
        for j in 0..(mesh_ky + 1){
            let kk = i_j_to_kk(i, j, mesh_kx, mesh_ky, false, system.size());

            let seud = diag(kk + deltakk,system).sort();

            seud_mat[i][j] = Some(seud);
        }
    }

    seud_mat
}
