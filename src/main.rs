use uuuddd2::{
    cal_bc_bcd::{get_tanzakus_rayon,}, calbinfo::cal_bcd_grid, model::{Param, System}
};

fn main() {
    let _comment = "a2";

    let system = System::UuudddTmd(Param::test());

    // let binfos_merged_onkks = calculate_band_info_all_band(false, 100, &system, comment);

    // let fermi_energy = -0.15;
    // println!(
    //     "bc_normal : {}\nbc_grid   : {}",
    //     bcd_sum_up(fermi_energy, &binfos_merged_onkks),
    //     cal_bcd_grid_test(fermi_energy, &system)
    // );

    let tanzaku = get_tanzakus_rayon(&system, 200, 2000, 0.3);

    tanzaku.export_integrated_bc_bcd("test.dat").expect("wow");

    println!("{}",cal_bcd_grid(0.0, &system, 10, 100, 0.1))
}
