use uuuddd2::{
    cal_bc_bcd::{get_tanzakus_rayon,MeshInfo}, calbinfo::calculate_band_info_all_band, model::{Param, System}
};

use std::time::Instant;

fn main() {
    let _comment = "a2";

    let start = Instant::now();

    let system = System::UuudddTmd(Param::sj());

    // let _binfos_merged_onkks = calculate_band_info_all_band(true, 2000, &system, comment);

    // let fermi_energy = -0.15;
    // println!(
    //     "bc_normal : {}\nbc_grid   : {}",
    //     bcd_sum_up(fermi_energy, &binfos_merged_onkks),
    //     cal_bcd_grid_test(fermi_energy, &system)
    // );

    let mesh_info = MeshInfo::new(400, 300, 0.2);

    let tanzaku = get_tanzakus_rayon(&system, &mesh_info,-12.1,12.1,true);

    tanzaku.export_integrated_bc_bcd(&system,&mesh_info).expect("except");

    let duration = start.elapsed();

    println!("実行時間: {:?}", duration)
}