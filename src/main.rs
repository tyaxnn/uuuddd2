use uuuddd2::{
    export::{calculate_band_info_all_band,},
    model::{System,Param},
};

fn main() {
    let comment = "a2";

    calculate_band_info_all_band(500, &System::UuudddTmd(Param::test()), comment);
}
