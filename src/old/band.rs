use crate::consts::*;
use crate::diag::{diag_uuuddd2,};
use crate::dv2::DV2;
use crate::berry::point_in_triangle_simple;

use nalgebra::{Vector2,};
use std::fs;
use std::fmt::Write;

//--------------------------------------------------------------------
//                          バンド図<->波数空間                    
//--------------------------------------------------------------------

//バンド図における横軸から対応するk空間のベクトルを返す
fn band_2_kk (x : f64) -> Vector2<f64> {
    if x < G_X_B{
        {
            let div = (x - MPX_B) / (G_X_B - MPX_B);
            MP_KS * (1. - div) + GAMMA * div
        }
    }
    else if  x < M_X_B{
        {
            let div = (x - G_X_B) / (M_X_B - G_X_B);
            GAMMA * (1. - div) + MINKS * div
        }
    }
    else if x < K_X_B{
        {
            let div = (x - M_X_B) / (K_X_B - M_X_B);
            MINKS * (1. - div) + KINKS * div
        }
    }
    else{
        {
            let div = (x - K_X_B) / (G2X_B - K_X_B);
            KINKS * (1. - div) + GAMMA * div
        }
    }
}

//--------------------------------------------------------------------
//                          処理                   
//--------------------------------------------------------------------

pub fn export_band_dat_uuuddd(lambda : f64, jj : f64, mu : f64, comment : &str) {
    let graph_div = 1000;

    let mut file_str = "x u1 u2 u3 u4 u5 u6 d1 d2 d3 d4 d5 d6".to_string();

    for i in 0..(graph_div + 1){
        let per = i as f64 / graph_div as f64;

        let x = per * MPX_B + (1. - per) * G2X_B;

        let seud = diag_uuuddd2(band_2_kk(x), lambda, jj, mu);

        file_str = format!(
            "{}\n{} {} {} {} {} {} {} {} {} {} {} {} {}",
            file_str,
            x,
            seud.u.eigenvalues[0],
            seud.u.eigenvalues[1],
            seud.u.eigenvalues[2],
            seud.u.eigenvalues[3],
            seud.u.eigenvalues[4],
            seud.u.eigenvalues[5],
            seud.d.eigenvalues[0],
            seud.d.eigenvalues[1],
            seud.d.eigenvalues[2],
            seud.d.eigenvalues[3],
            seud.d.eigenvalues[4],
            seud.d.eigenvalues[5],
        );
    }
    fs::write(format!("./output/bands/dats/kk/band_kk_lambda{}_j{}_mu{}_{}.dat",lambda,jj,mu,comment), file_str).unwrap();

}

pub fn export_band_dat_uuuddd_2d(graph_mesh : usize,lambda : f64, jj : f64, mu : f64, comment : &str){
    let mesh_kx = graph_mesh;
    let mesh_ky  = graph_mesh;

    let dv2_k1 = DV2::from_car(KPPKS) - DV2::from_car(-KINKS);
    let dv2_k2 = DV2::from_car(KP_KS) - DV2::from_car(-KINKS);

    let mut file_str = format!("kx ky u1 u2 u3 u4 u5 u6 d1 d2 d3 d4 d5 d6 \ngraphmesh {} lambda {} jj {} mu {}",graph_mesh,lambda, jj, mu);

    for i in 0..(mesh_kx){
        for j in 0..(mesh_ky){
            let if64 = i as f64 / mesh_kx as f64;
            let jf64 = j as f64 / mesh_ky as f64;

            let mut kk_dv2 = dv2_k1 * if64 + dv2_k2 * jf64 + DV2::from_car(-KINKS);

            let seud = diag_uuuddd2(kk_dv2.to_car(),lambda,jj,mu).sort();

            if point_in_triangle_simple(kk_dv2.to_car(), KINKS, 2. * KINKS, KPPKS){
                kk_dv2 = kk_dv2 - dv2_k1;
            }
            else if point_in_triangle_simple(kk_dv2.to_car(), KINKS, 2. * KINKS, KP_KS){
                kk_dv2 = kk_dv2 - dv2_k2;
            }

            let kx = kk_dv2.to_car().x;
            let ky = kk_dv2.to_car().y;

            let new_line = format!(
                "{} {} {} {} {} {} {} {} {} {} {} {} {} {}",
                kx,ky,
                seud.u.eigenvalues[0],
                seud.u.eigenvalues[1],
                seud.u.eigenvalues[2],
                seud.u.eigenvalues[3],
                seud.u.eigenvalues[4],
                seud.u.eigenvalues[5],
                seud.d.eigenvalues[0],
                seud.d.eigenvalues[1],
                seud.d.eigenvalues[2],
                seud.d.eigenvalues[3],
                seud.d.eigenvalues[4],
                seud.d.eigenvalues[5],
            );

            write!(file_str, "\n{}", new_line).unwrap();

        }
    }
    fs::write(format!("./output/bands/dats/2d/band_2d_lambda{}_j{}_mu{}_{}.dat",lambda,jj,mu,comment), file_str).unwrap();

}