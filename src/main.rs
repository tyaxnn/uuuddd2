use uuuddd2::{
    calbinfo::{calculate_band_info_all_band,cal_bc_grid_test,cal_bcd_grid_test},
    cal_bc_bcd::{bc_sum_up,bcd_sum_up,Tanzakus,get_tanzakus},
    model::{System,Param},
};

fn main() {
    let comment = "a2";

    let system = System::UuudddTmd(Param::test());

    // let binfos_merged_onkks = calculate_band_info_all_band(false, 100, &system, comment);

    // let fermi_energy = -0.15;
    // println!(
    //     "bc_normal : {}\nbc_grid   : {}",
    //     bcd_sum_up(fermi_energy, &binfos_merged_onkks),
    //     cal_bcd_grid_test(fermi_energy, &system)
    // );

    let tanzaku = get_tanzakus(&system, 10, 100, 0.1);

    tanzaku.export_integrated_bc_bcd("test.dat").expect("wow");
}
