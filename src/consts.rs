use nalgebra::{Complex,Vector2,};

//--------------------------------------------------------------------
//                          定数の定義                    
//--------------------------------------------------------------------

pub const SQRT_3 : f64 = 1.732050807568877293527446341505872367_f64;
pub const PI     : f64 = std::f64::consts::PI;

pub const I      : Complex<f64> = Complex::new(0.,1.);
pub const ZERO   : Complex<f64> = Complex::new(0.,0.);
pub const ONE    : Complex<f64> = Complex::new(1.,0.);

//--------------------------------------------------------------------
//                          実逆格子のベクトルを定義                    
//--------------------------------------------------------------------

//最近接サイト間の長さ
const A : f64 = 1.0;

//格子空間の次近接ベクトル
pub const A1 :  Vector2<f64> = Vector2::new(A * -0.5 * SQRT_3, A * SQRT_3 / 2. * SQRT_3);
pub const A2 :  Vector2<f64> = Vector2::new(A * -0.5 * SQRT_3, A * -1. * SQRT_3 / 2. * SQRT_3);
pub const A3 :  Vector2<f64> = Vector2::new(A * SQRT_3,0.);

//格子空間の最近接ベクトル
pub const D1 :  Vector2<f64> = Vector2::new(A * SQRT_3 / 2.,A * 0.5);
pub const D2 :  Vector2<f64> = Vector2::new(A * -1. * SQRT_3 / 2.,A * 0.5);
pub const D3 :  Vector2<f64> = Vector2::new(0.,A * -1.);

const TRI : f64 = 4. * PI / 9. / A;

//k空間における特徴点
pub const KP_KS : Vector2<f64> = Vector2::new(TRI * SQRT_3 / 2.,TRI * 0.5 );
pub const KPPKS : Vector2<f64> = Vector2::new(0.,-TRI);
pub const GAMMA : Vector2<f64> = Vector2::new(0.,0.);
pub const KINKS : Vector2<f64> = Vector2::new(TRI * SQRT_3 / 2.,-TRI * 0.5 );
pub const MINKS : Vector2<f64> = Vector2::new(TRI * SQRT_3 / 2.,0. );
pub const MP_KS : Vector2<f64> = Vector2::new(-TRI * SQRT_3 / 2.,0. );

pub const KP_KS2 : Vector2<f64> = Vector2::new(4. * SQRT_3 * PI / 9. * 0.5 / A,4. * SQRT_3 * PI / 9. * SQRT_3 * 0.5 / A);
pub const KPPKS2 : Vector2<f64> = Vector2::new(4. * SQRT_3 * PI / 9. * 0.5 / A,-4. * SQRT_3 * PI / 9. * SQRT_3 * 0.5 / A);
pub const GAMMA2 : Vector2<f64> = Vector2::new(0.,0.);
pub const KINKS2 : Vector2<f64> = Vector2::new(4. * SQRT_3 * PI / 9. / A,0.);
pub const MINKS2 : Vector2<f64> = Vector2::new(SQRT_3 * PI / 3. / A,PI / 3. / A);

//バンド図におけるx軸の特徴点
pub const MPX_B : f64 = -TRI * SQRT_3 / 2.;
pub const G_X_B : f64 = 0.;
pub const M_X_B : f64 = TRI * SQRT_3 / 2.;
pub const K_X_B : f64 = (0.5 + SQRT_3/2.) * TRI;
pub const G2X_B : f64 = (1.5 + SQRT_3/2.) * TRI;

//--------------------------------------------------------------------
//                          ハミルトニアンのパラメーター                    
//--------------------------------------------------------------------

//hopping
pub const T      : f64 = 1.0;