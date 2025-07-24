use nalgebra::Vector2;

use crate::binfo::BinfosMergedOnkk;
use crate::util::GridInfo;
use crate::model::System;
use crate::calbinfo::calculate_band_info_grid;

use std::fs::File;
use std::io::{Write, BufWriter};
use rayon::prelude::*;

pub struct Tanzakus{
    data : Vec<Tanzaku>,
    ground_energy : f64,
    max_energy : f64,
    div_tanzaku : usize,
}

impl Tanzakus{
    fn new(ground_energy : f64, max_energy : f64, div_tanzaku : usize) -> Self{
        let data = Tanzaku::ini_vec(ground_energy, max_energy, div_tanzaku);

        Tanzakus { data, ground_energy, max_energy, div_tanzaku}
    }
    pub fn export_integrated_bc_bcd(&self, filename: &str) -> std::io::Result<()> {
        let file = File::create(format!("./output/tanzaku/{}",filename))?;
        let mut writer = BufWriter::new(file);

        writeln!(writer, "# energy,bc_sum,bcd_x_sum,bcd_y_sum")?;

        for i in 0..=self.div_tanzaku {
            let energy = self.ground_energy
                + (self.max_energy - self.ground_energy) * (i as f64) / (self.div_tanzaku as f64);

            let mut bc_sum = 0.0;
            let mut bcd_sum = Vector2::new(0.0, 0.0);

            for tanzaku in &self.data {
                // エネルギーが E 以下にある tanzaku を合計
                if tanzaku.h_energy <= energy {
                    bc_sum += tanzaku.bc;
                    bcd_sum += tanzaku.bcd;
                }
            }

            writeln!(
                writer,
                "{:.6},{:.6},{:.6},{:.6}",
                energy,
                bc_sum,
                bcd_sum.x,
                bcd_sum.y
            )?;
        }

        Ok(())
    }
    pub fn merge(&mut self, other: Tanzakus) {
        self.data.extend(other.data);
    }
}

pub struct Tanzaku{
    h_energy : f64,
    bc       : f64,
    bcd      : Vector2<f64>
}
impl Tanzaku{
    fn ini_vec(ground_energy : f64, max_energy : f64, div_tanzaku : usize) -> Vec<Tanzaku>{

        let mut out = Vec::new();
        for i in 0..div_tanzaku{
            let ip1_0_to_1 = i as f64 / div_tanzaku as f64;

            let h_energy = ground_energy + (max_energy - ground_energy) * ip1_0_to_1;

            let tanzaku = Tanzaku{h_energy,bc : 0.0, bcd : Vector2::zeros()};

            out.push(tanzaku)
 
        }

        out
    }
}


pub fn bc_tanzaku(tanzakus : &mut Tanzakus, binfos_merged_onkks : &Vec<BinfosMergedOnkk>){

    let tanzaku_energy = (tanzakus.max_energy - tanzakus.ground_energy) / (tanzakus.div_tanzaku as f64);

    for binfos_merged_onkk in binfos_merged_onkks{
        
        for binfo in &binfos_merged_onkk.infos{
            let ei_num = (((binfo.eigen - tanzakus.ground_energy))/ tanzaku_energy) as usize;

            tanzakus.data[ei_num].bc += binfo.berry;
            tanzakus.data[ei_num].bcd = tanzakus.data[ei_num].bcd + binfo.bcd.unwrap();
        }
    }
}

pub fn get_tanzakus(system : &System,grid_mesh : usize, max_sub_mesh : usize, mesh_scale : f64) -> Tanzakus{

    let mut tanzakus = Tanzakus::new(-3.0,3.10,1000);

    let n = grid_mesh;

    for grid_x in 0..n{
        for grid_y in 0..n{
            let grid_info = GridInfo::new_ijn(grid_x,grid_y,n,n);

            let (binfos_merged_onkks, max_bcd ) = calculate_band_info_grid(1, system, grid_info);

            if max_bcd.abs() > 1. / mesh_scale{

                let graph_mesh = std::cmp::min((max_bcd.abs() * mesh_scale) as usize,max_sub_mesh);
                let (binfos_merged_onkks_new, _ ) = calculate_band_info_grid(graph_mesh, system, grid_info);

                bc_tanzaku(&mut tanzakus, &binfos_merged_onkks_new);
            }
            else{
                bc_tanzaku(&mut tanzakus, &binfos_merged_onkks);
            }
        }
    }

    tanzakus
}

pub fn get_tanzakus_rayon(system : &System,grid_mesh : usize, max_sub_mesh : usize, mesh_scale : f64) -> Tanzakus{

    let n = grid_mesh;
    let total_points = n * n;

    // 全インデックスを準備
    let all_indices: Vec<(usize, usize)> = (0..n).flat_map(|i| {
        (0..n).map(move |j| (i, j))
    }).collect();

    // スレッド数程度にchunk分割（例：16スレッドなら16分割）
    let chunk_size = total_points / 100;

    let sub_tanzakus: Vec<Tanzakus> = all_indices
        .par_chunks(chunk_size)
        .map(|chunk| {
            let mut local_tanzakus = Tanzakus::new(-3.0, 3.10, 1000);
            for &(grid_x, grid_y) in chunk {
                let grid_info = GridInfo::new_ijn(grid_x, grid_y, n, n);
                let (binfos_merged_onkks, max_bcd) = calculate_band_info_grid(1, system, grid_info);

                let binfos = if max_bcd.abs() > 1. / mesh_scale {
                    let graph_mesh = std::cmp::min((max_bcd.abs() * mesh_scale) as usize, max_sub_mesh);
                    let (refined_binfos, _) = calculate_band_info_grid(graph_mesh, system, grid_info);
                    refined_binfos
                } else {
                    binfos_merged_onkks
                };

                bc_tanzaku(&mut local_tanzakus, &binfos);
            }
            local_tanzakus
        })
        .collect();

    // 最終統合（スレッド数回だけmerge）
    let mut final_tanzakus = Tanzakus::new(-3.0, 3.10, 1000);
    for t in sub_tanzakus {
        final_tanzakus.merge(t);
    }

    final_tanzakus
}

pub fn bc_sum_up(fermi_energy : f64 ,binfos_merged_onkks : &Vec<BinfosMergedOnkk>) -> f64 {

    let mut bc = 0.;
    
    for binfos_merged_onkk in binfos_merged_onkks{
        
        for binfo in &binfos_merged_onkk.infos{
            if binfo.eigen < fermi_energy{
                bc += binfo.berry;
            }
        }
    }

    bc 
}

pub fn bcd_sum_up(fermi_energy : f64 ,binfos_merged_onkks : &Vec<BinfosMergedOnkk>) -> Vector2<f64> {

    let mut bcd = Vector2::zeros();
    
    for binfos_merged_onkk in binfos_merged_onkks{
        
        for binfo in &binfos_merged_onkk.infos{
            if binfo.eigen < fermi_energy{
                bcd = bcd + binfo.bcd.unwrap();
            }
        }
    }

    bcd
}