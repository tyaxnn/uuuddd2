use crate::consts::T;

pub enum System{
    Uuuddd(Param),
    UuudddTmd(Param),
    Sato(Param),
    SatoTmd(Param),
    Tmd(Param),
    UuudddTmdUM(Param),
}

impl System{
    pub fn size(&self) -> usize{
        match self{
            Self::Uuuddd(_) => {6},
            Self::UuudddTmd(_) => {6}
            Self::Sato(_) => {2}
            Self::SatoTmd(_) => {2}
            Self::Tmd(_) => {2}
            Self::UuudddTmdUM(_) => {6}
        }
    }
    pub fn debug(&self) -> String{
        match self {
            Self::Uuuddd(param) => {
                format!("uuuddd_lambda{}_j{}_mu{}",
                    param.lambda.to_string().replace('.', "p")
                    ,param.jj.to_string().replace('.', "p")
                    ,param.mu.to_string().replace('.', "p")
                )
            },
            Self::UuudddTmd(param) => {
                format!("tmduuuddd_lambda{}_j{}_mu{}",
                    param.lambda.to_string().replace('.', "p")
                    ,param.jj.to_string().replace('.', "p")
                    ,param.mu.to_string().replace('.', "p")
                )
            }
            Self::Sato(param) => {
                format!("sato_lambda{}_j{}_mu{}",
                    param.lambda.to_string().replace('.', "p")
                    ,param.jj.to_string().replace('.', "p")
                    ,param.mu.to_string().replace('.', "p")
                )
            }
            Self::SatoTmd(param) => {
                format!("tmdsato_lambda{}_j{}_mu{}",
                    param.lambda.to_string().replace('.', "p")
                    ,param.jj.to_string().replace('.', "p")
                    ,param.mu.to_string().replace('.', "p")
                )
            }
            Self::Tmd(param) => {
                format!("tmd_lambda{}_j{}_mu{}",
                    param.lambda.to_string().replace('.', "p")
                    ,param.jj.to_string().replace('.', "p")
                    ,param.mu.to_string().replace('.', "p")
                )
            }
            Self::UuudddTmdUM(param) => {
                format!("tmduuudddum_lambda{}_j{}_mu{}",
                    param.lambda.to_string().replace('.', "p")
                    ,param.jj.to_string().replace('.', "p")
                    ,param.mu.to_string().replace('.', "p")
                )
            }
        }
    }
}

pub struct Param{
    pub lambda : f64,
    pub jj : f64,
    pub mu : f64,
}

impl Param{
    pub fn new(lambda: f64, jj: f64, mu: f64) -> Self{
        Param { lambda, jj, mu }
    }
    pub fn test() -> Self{
        let lambda = 0.1 * T;
        let jj = 0.25;
        let mu = 0.0;

        Param { lambda, jj, mu }
    }
    pub fn sj() -> Self{
        let lambda = 0.1 * T;
        let jj = 10.0;
        let mu = 0.0;

        Param { lambda, jj, mu }
    }
}