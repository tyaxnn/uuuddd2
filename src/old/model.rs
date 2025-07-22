pub struct Model{
    pub lambda : f64,       //次近接hopping係数
    pub jj : f64,           //spin秩序の大きさ
    pub mu : f64,           //化学ポテンシャル
    pub tmd : f64,          //次近接ホッピングの異方性
    pub hopping : Hopping,  //Kane-Mele or Modified-Kane-Mele
    pub order : Order       //spin秩序
}

impl Model{
    pub fn uuuddd(lambda: f64, jj: f64, mu: f64) -> Self{
        Model { lambda , jj , mu , tmd: 1., hopping: Hopping::Mod, order : Order::from_str("uuuddd").unwrap() }
    }
}

pub enum Hopping{
    Kan,
    Mod,
}

pub struct Order{
    a : Spin,
    b : Spin,
    c : Spin,
    d : Spin,
    e : Spin,
    f : Spin,
}

impl Order {
    pub fn from_str(input: &str) -> Result<Self, String> {
        if input.len() != 6 {
            return Err("入力文字列は6文字である必要があります。".to_string());
        }

        let chars: Vec<char> = input.chars().collect();

        let to_spin = |c: char| -> Result<Spin, String> {
            match c {
                'u' => Ok(Spin::u),
                'd' => Ok(Spin::d),
                _ => Err(format!("無効な文字 '{}' が含まれています。", c)),
            }
        };

        Ok(Order {
            a: to_spin(chars[0])?,
            b: to_spin(chars[1])?,
            c: to_spin(chars[2])?,
            d: to_spin(chars[3])?,
            e: to_spin(chars[4])?,
            f: to_spin(chars[5])?,
        })
    }
}

pub enum Spin{
    u,
    d,
}