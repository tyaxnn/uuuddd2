use crate::consts::*;
use crate::diag::{diag,SEudEnum};
use crate::util::{i_j_to_kk,cal_cell_area};
use crate::model::System;
use crate::binfo::BinfosMergedOnkk;

use nalgebra::{Vector2,};
use std::fs;

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

    let seud_mat = make_seud_mat(Vector2::zeros(), mesh_kx, mesh_ky, system);
    let seud_mat_px = make_seud_mat(Vector2::new(DELTA,0.0), mesh_kx, mesh_ky, system);
    let seud_mat_py = make_seud_mat(Vector2::new(0.0,DELTA), mesh_kx, mesh_ky, system);

    let cell_area = cal_cell_area(mesh_kx, mesh_ky, system.size());

    for i in 0..(mesh_kx){
        for j in 0..(mesh_ky){
            let kk = i_j_to_kk(i, j, mesh_kx, mesh_ky, true, system.size());
            
            let binfos_merged_onkk = BinfosMergedOnkk::cal(i,j,&seud_mat,&seud_mat_px,&seud_mat_py,system.size(),cell_area,kk);

            binfos_merged_onkk.write(&mut file_str, system.size());
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
