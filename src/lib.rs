pub mod consts;         //定数を定義する
mod diag;               //システムごとの対角化関数
mod dv2;            //福井-初貝のgridを菱形に取る際の菱形を分かりやすく示す構造体を定義する
pub mod calbinfo;       //k点ごとの全てのバンドの計算したい量を一括で計算する関数を定義する
mod util;               //細々とした関数の詰め合わせ
pub mod model;          //システムを定義する
mod binfo;              //バンドの情報を格納する構造体を定義する
pub mod cal_bc_bcd;