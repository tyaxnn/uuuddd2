use crate::consts::*;
use crate::dv2::DV2;

use nalgebra::{Vector2,};

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


pub fn i_j_to_kk(i : usize, j : usize, mesh_kx : usize, mesh_ky : usize, hex : bool, size : usize) -> Vector2<f64>{
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

            let kx = kk_dv2.to_car().x;
            let ky = kk_dv2.to_car().y;

            Vector2::new(kx,ky)
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

            let kx = kk_dv2.to_car2().x;
            let ky = kk_dv2.to_car2().y;

            Vector2::new(kx,ky)
        }
        _ => {panic!("invalid size")}
    }
}

pub fn cal_cell_area(mesh_kx : usize, mesh_ky : usize, size : usize) -> f64{
    match size {
        2 => {
            let dv2_k1 = DV2::from_car2(KPPKS2) - DV2::from_car2(-KINKS2);
            let dv2_k2 = DV2::from_car2(KP_KS2) - DV2::from_car2(-KINKS2);

            let lattice_1_len = dv2_k1.to_car2().norm() / mesh_kx as f64;
            let lattice_2_len = dv2_k2.to_car2().norm() / mesh_ky as f64;

            let cell_area = lattice_1_len * lattice_2_len * SQRT_3 * 0.5;

            cell_area
        }
        6 => {
            let dv2_k1 = DV2::from_car(KPPKS) - DV2::from_car(-KINKS);
            let dv2_k2 = DV2::from_car(KP_KS) - DV2::from_car(-KINKS);

            let lattice_1_len = dv2_k1.to_car().norm() / mesh_kx as f64;
            let lattice_2_len = dv2_k2.to_car().norm() / mesh_ky as f64;

            let cell_area = lattice_1_len * lattice_2_len * SQRT_3 * 0.5;

            cell_area
        }
        _ => panic!("invalid size")
    }
}