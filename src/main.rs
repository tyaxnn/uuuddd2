use uuuddd2::{
    export::{calculate_band_info_all_band,},
    model::{System,Param},
};

fn main() {
    let comment = "a1";

    calculate_band_info_all_band(true, 500, &System::UuudddTmdUM(Param::test()), comment);
}
