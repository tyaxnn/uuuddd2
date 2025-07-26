use crate::consts::*;
use crate::diag::{SEudEnum,};

use nalgebra::{Complex, Vector2, Vector6,};
use std::fmt::Write;

//波数点とその上の全てのバンドにおけるBinfoを格納する構造体
pub struct BinfosMergedOnkk{
    kk : Vector2<f64>,
    pub infos : Vec<Binfo>
}

impl BinfosMergedOnkk{
    pub fn cal_cd(
        i : usize, 
        j : usize, 
        seud_mat : &Vec<Vec<Option<SEudEnum>>>,
        seud_mat_px : &Vec<Vec<Option<SEudEnum>>>,
        seud_mat_py : &Vec<Vec<Option<SEudEnum>>>,
        seud_mat_mx : &Vec<Vec<Option<SEudEnum>>>,
        seud_mat_my : &Vec<Vec<Option<SEudEnum>>>,
        size : usize,
        kk : Vector2<f64>
    ) -> Self{
        let binfos_sepud = BinfosSEPud::cal_cd(i, j, seud_mat, seud_mat_px, seud_mat_py,seud_mat_mx, seud_mat_my, size);

        let binfos_merged_onkk = binfos_sepud.mergedsorted(kk);

        binfos_merged_onkk        
    }
    //構造体を計算して作る関数
    pub fn cal(
        i : usize, 
        j : usize, 
        seud_mat : &Vec<Vec<Option<SEudEnum>>>,
        seud_mat_px : &Vec<Vec<Option<SEudEnum>>>,
        seud_mat_py : &Vec<Vec<Option<SEudEnum>>>,
        size : usize,
        kk : Vector2<f64>
    ) -> Self{
        let binfos_sepud = BinfosSEPud::cal(i, j, seud_mat, seud_mat_px, seud_mat_py, size);

        let binfos_merged_onkk = binfos_sepud.mergedsorted(kk);

        binfos_merged_onkk        
    }
    //構造体の情報を入力のStringの最後の行にS加える関数
    pub fn write(&self, file_str : &mut String, size : usize, cell_area : f64){

        let kk = self.kk;
        let berrys : Vec<f64> = self.infos.iter().map(|b| b.berry).collect();
        let eigens : Vec<f64> = self.infos.iter().map(|b| b.eigen).collect();
        let spins  : Vec<u8> = self.infos.iter().map(|b| b.spin.to_u8()).collect();
        let bcds   : Vec<Vector2<f64>> = self.infos.iter().map(|b| b.bcd.unwrap()).collect();

        write!(file_str,"{},{}",kk.x,kk.y).unwrap();

        for k in 0..size{
            write!(file_str,",{}",
                berrys[k] / cell_area,
            ).unwrap();
        }
        for k in 0..size{
            write!(file_str,",{}",
                eigens[k],
            ).unwrap();
        }
        for k in 0..size{
            write!(file_str,",{}",
                spins[k],
            ).unwrap();
        }
        for k in 0..size{
            write!(file_str,",{},{}",
                bcds[k].x,bcds[k].y
            ).unwrap();
        }

        write!(file_str,"\n").unwrap();
    }
    pub fn max_bc(&self, cell_area : f64) -> f64 {
        
        let berrys : Vec<f64> = self.infos.iter().map(|b| b.berry).collect();

        let max_berry = berrys
            .iter()
            .copied()        // Option<&f64> → Option<f64>
            .max_by(
                |a, b| 
                a.abs().partial_cmp(&b.abs()).unwrap()
            );

        max_berry.unwrap() / cell_area
    }
    pub fn max_bcd(&self, cell_area : f64) -> f64 {
        
        let bcds : Vec<f64> = self.infos.iter().map(|b| b.bcd.unwrap().norm()).collect();

        let max_bcd = bcds
            .iter()
            .copied()        // Option<&f64> → Option<f64>
            .max_by(
                |a, b| 
                a.abs().partial_cmp(&b.abs()).unwrap()
            );

        max_bcd.unwrap() / cell_area
    }
}

struct BinfosSEPud{
    u : Vec<Binfo>,
    d : Vec<Binfo>
}
impl BinfosSEPud {
    pub fn mergedsorted(&self,kk : Vector2<f64>) -> BinfosMergedOnkk {
        let mut infos = [self.u.clone(), self.d.clone()].concat();

        infos.sort_by(|a, b| a.eigen.partial_cmp(&b.eigen).unwrap());

        BinfosMergedOnkk { kk, infos }
    }
    fn cal_cd(
        i : usize, 
        j : usize, 
        seud_mat : &Vec<Vec<Option<SEudEnum>>>,
        seud_mat_px : &Vec<Vec<Option<SEudEnum>>>,
        seud_mat_py : &Vec<Vec<Option<SEudEnum>>>,
        seud_mat_mx : &Vec<Vec<Option<SEudEnum>>>,
        seud_mat_my : &Vec<Vec<Option<SEudEnum>>>,
        size : usize,
    ) -> Self{
        let mut binfos_sepud = BinfosSEPud::cal_without_bcd(i, j, seud_mat,size);

        let binfos_sepud_px = BinfosSEPud::cal_without_bcd(i, j, seud_mat_px,size);
        let binfos_sepud_py = BinfosSEPud::cal_without_bcd(i, j, seud_mat_py,size);
        let binfos_sepud_mx = BinfosSEPud::cal_without_bcd(i, j, seud_mat_mx,size);
        let binfos_sepud_my = BinfosSEPud::cal_without_bcd(i, j, seud_mat_my,size);

        for ei in 0..size{
            let berry_px = binfos_sepud_px.u[ei].berry;
            let berry_py = binfos_sepud_py.u[ei].berry;
            let berry_mx = binfos_sepud_mx.u[ei].berry;
            let berry_my = binfos_sepud_my.u[ei].berry;

            binfos_sepud.u[ei].bcd = Some({
                let changex = (berry_px - berry_mx) / DELTA / 2.;
                let changey = (berry_py - berry_my) / DELTA / 2.;

                Vector2::new(changex,changey)
            });
        }
        for ei in 0..size{
            let berry_px = binfos_sepud_px.d[ei].berry;
            let berry_py = binfos_sepud_py.d[ei].berry;
            let berry_mx = binfos_sepud_mx.d[ei].berry;
            let berry_my = binfos_sepud_my.d[ei].berry;

            binfos_sepud.d[ei].bcd = Some({
                let changex = (berry_px - berry_mx) / DELTA / 2.;
                let changey = (berry_py - berry_my) / DELTA / 2.;

                Vector2::new(changex,changey)
            });
        }   

        binfos_sepud
    }
    fn cal(
        i : usize, 
        j : usize, 
        seud_mat : &Vec<Vec<Option<SEudEnum>>>,
        seud_mat_px : &Vec<Vec<Option<SEudEnum>>>,
        seud_mat_py : &Vec<Vec<Option<SEudEnum>>>,
        size : usize,
    ) -> Self{
        let mut binfos_sepud = BinfosSEPud::cal_without_bcd(i, j, seud_mat,size);

        let binfos_sepud_px = BinfosSEPud::cal_without_bcd(i, j, seud_mat_px,size);
        let binfos_sepud_py = BinfosSEPud::cal_without_bcd(i, j, seud_mat_py,size);

        for ei in 0..size{
            let berry = binfos_sepud.u[ei].berry;
            let berry_px = binfos_sepud_px.u[ei].berry;
            let berry_py = binfos_sepud_py.u[ei].berry;

            binfos_sepud.u[ei].bcd = Some({
                let changex = (berry_px - berry) / DELTA;
                let changey = (berry_py - berry) / DELTA;

                Vector2::new(changex,changey)
            });
        }
        for ei in 0..size{
            let berry = binfos_sepud.d[ei].berry;
            let berry_px = binfos_sepud_px.d[ei].berry;
            let berry_py = binfos_sepud_py.d[ei].berry;

            binfos_sepud.d[ei].bcd = Some({
                let changex = (berry_px - berry) / DELTA;
                let changey = (berry_py - berry) / DELTA;

                Vector2::new(changex,changey)
            });
        }   

        binfos_sepud
    }
    fn cal_without_bcd(
        i : usize, 
        j : usize, 
        seud_mat : &Vec<Vec<Option<SEudEnum>>>,
        size : usize,  
    ) -> Self{
        match size{
            2 => {

                let seud = SEudEnum::get_cont(&seud_mat[i][j]).is_2();
                let seud_pi = SEudEnum::get_cont(&seud_mat[i+1][j]).is_2();
                let seud_pij = SEudEnum::get_cont(&seud_mat[i+1][j+1]).is_2();
                let seud_pj = SEudEnum::get_cont(&seud_mat[i][j+1]).is_2();

                let mut berrys_up = vec![Binfo::ini();2];
                let mut berrys_do = vec![Binfo::ini();2];

                for ei in 0..2{
                    let u00 : Vector2<Complex<f64>> = seud.u.eigenvectors.column(ei).into();
                    let u10 : Vector2<Complex<f64>> = seud_pi.u.eigenvectors.column(ei).into();
                    let u11 : Vector2<Complex<f64>> = seud_pij.u.eigenvectors.column(ei).into();
                    let u01 : Vector2<Complex<f64>> = seud_pj.u.eigenvectors.column(ei).into(); 

                    let berry = {
                        let u1 = u00.dotc(&u10);
                        let u2 = u10.dotc(&u11);
                        let u3 = u11.dotc(&u01);
                        let u4 = u01.dotc(&u00);

                        (u1 * u2 * u3 * u4).arg()
                    } ;

                    berrys_up[ei] = Binfo::new(berry, seud.u.eigenvalues[ei],Spin::U,None)

                }
                for ei in 0..2{
                    let u00 : Vector2<Complex<f64>> = seud.d.eigenvectors.column(ei).into();
                    let u10 : Vector2<Complex<f64>> = seud_pi.d.eigenvectors.column(ei).into();
                    let u11 : Vector2<Complex<f64>> = seud_pij.d.eigenvectors.column(ei).into();
                    let u01 : Vector2<Complex<f64>> = seud_pj.d.eigenvectors.column(ei).into(); 

                    let berry = {
                        let u1 = u00.dotc(&u10);
                        let u2 = u10.dotc(&u11);
                        let u3 = u11.dotc(&u01);
                        let u4 = u01.dotc(&u00);

                        (u1 * u2 * u3 * u4).arg()
                    } ;

                    berrys_do[ei] = Binfo::new(berry, seud.d.eigenvalues[ei],Spin::D,None)

                }

                BinfosSEPud{u : berrys_up, d : berrys_do}
            }
            6 => {
                let seud = SEudEnum::get_cont(&seud_mat[i][j]).is_6();
                let seud_pi = SEudEnum::get_cont(&seud_mat[i+1][j]).is_6();
                let seud_pij = SEudEnum::get_cont(&seud_mat[i+1][j+1]).is_6();
                let seud_pj = SEudEnum::get_cont(&seud_mat[i][j+1]).is_6();

                let mut berrys_up = vec![Binfo::ini();6];
                let mut berrys_do = vec![Binfo::ini();6];

                for ei in 0..6{
                    let u00 : Vector6<Complex<f64>> = seud.u.eigenvectors.column(ei).into();
                    let u10 : Vector6<Complex<f64>> = seud_pi.u.eigenvectors.column(ei).into();
                    let u11 : Vector6<Complex<f64>> = seud_pij.u.eigenvectors.column(ei).into();
                    let u01 : Vector6<Complex<f64>> = seud_pj.u.eigenvectors.column(ei).into(); 

                    let berry = {
                        let u1 = u00.dotc(&u10);
                        let u2 = u10.dotc(&u11);
                        let u3 = u11.dotc(&u01);
                        let u4 = u01.dotc(&u00);

                        (u1 * u2 * u3 * u4).arg()
                    } ;

                    berrys_up[ei] = Binfo::new(berry, seud.u.eigenvalues[ei],Spin::U,None)

                }
                for ei in 0..6{
                    let u00 : Vector6<Complex<f64>> = seud.d.eigenvectors.column(ei).into();
                    let u10 : Vector6<Complex<f64>> = seud_pi.d.eigenvectors.column(ei).into();
                    let u11 : Vector6<Complex<f64>> = seud_pij.d.eigenvectors.column(ei).into();
                    let u01 : Vector6<Complex<f64>> = seud_pj.d.eigenvectors.column(ei).into(); 

                    let berry = {
                        let u1 = u00.dotc(&u10);
                        let u2 = u10.dotc(&u11);
                        let u3 = u11.dotc(&u01);
                        let u4 = u01.dotc(&u00);

                        (u1 * u2 * u3 * u4).arg()
                    } ;

                    berrys_do[ei] = Binfo::new(berry, seud.d.eigenvalues[ei],Spin::D,None)

                }

                BinfosSEPud{u : berrys_up, d : berrys_do}
            }
            _ => panic!("invalid size")
        }        
    }
}


#[derive(Debug, Clone, Copy)]
pub struct Binfo{
    pub berry : f64,
    pub eigen : f64,
    pub spin : Spin,
    pub bcd : Option<Vector2<f64>>
}

impl Binfo{
    fn ini() -> Self{
        Binfo { berry: 0.0, eigen: 0.0, spin : Spin::Undefined, bcd : None }
    }
    fn new(berry : f64, eigen : f64, spin : Spin, bcd : Option<Vector2<f64>>) -> Self{
        Binfo { berry, eigen , spin, bcd}
    }
}

#[derive(Debug, Clone, Copy)]
pub enum Spin{
    U,
    D,
    Undefined,
}

impl Spin{
    fn to_u8(self) -> u8{
        match self{
            Self::U => {1}
            Self::D => {0}
            Self::Undefined => {100}
        }
    }
}