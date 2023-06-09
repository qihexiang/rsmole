use strum_macros::EnumString;

#[derive(EnumString, Clone, Copy, PartialEq, Eq, Debug)]
pub enum Element {
    H = 1,
    He,
    Li,
    Be,
    B,
    C,
    N,
    O,
    F,
    Ne,
    Na,
    Mg,
    Al,
    Si,
    P,
    S,
    Cl,
    Ar,
    K,
    Ca,
    Sc,
    Ti,
    V,
    Cr,
    Mn,
    Fe,
    Co,
    Ni,
    Cu,
    Zn,
    Ga,
    Ge,
    As,
    Se,
    Br,
    Kr,
    Rb,
    Sr,
    Y,
    Zr,
    Nb,
    Mo,
    Tc,
    Ru,
    Rh,
    Pd,
    Ag,
    Cd,
    In,
    Sn,
    Sb,
    Te,
    I,
    Xe,
    Cs,
    Ba,
    La,
    Ce,
    Pr,
    Nd,
    Pm,
    Sm,
    Eu,
    Gd,
    Tb,
    Dy,
    Ho,
    Er,
    Tm,
    Yb,
    Lu,
    Hf,
    Ta,
    W,
    Re,
    Os,
    Ir,
    Pt,
    Au,
    Hg,
    Tl,
    Pb,
    Bi,
    Po,
    At,
    Rn,
    Fr,
    Ra,
    Ac,
    Th,
    Pa,
    U,
    Np,
    Pu,
    Am,
    Cm,
    Bk,
    Cf,
    Es,
    Fm,
    Md,
    No,
    Lr,
    Rf,
    Db,
    Sg,
    Bh,
    Hs,
    Mt,
    Ds,
    Rg,
    Cn,
    Nh,
    Fl,
    Mc,
    Lv,
    Ts,
    Og,
}

impl Element {
    pub fn get_average_mass(&self) -> f64 {
        match self {
            Element::H => 1.0080,
            Element::He => 4.0026,
            Element::Li => 6.94,
            Element::Be => 9.0122,
            Element::B => 10.81,
            Element::C => 12.011,
            Element::N => 14.007,
            Element::O => 15.999,
            Element::F => 18.998,
            Element::Ne => 20.180,
            Element::Na => 22.990,
            Element::Mg => 24.305,
            Element::Al => 26.982,
            Element::Si => 28.0855,
            Element::P => 30.974,
            Element::S => 32.06,
            Element::Cl => 35.45,
            Element::Ar => 39.95,
            Element::K => 39.098,
            Element::Ca => 40.078,
            Element::Sc => 44.956,
            Element::Ti => 47.867,
            Element::V => 50.942,
            Element::Cr => 51.996,
            Element::Mn => 54.938,
            Element::Fe => 55.845,
            Element::Co => 58.933,
            Element::Ni => 58.693,
            Element::Cu => 63.546,
            Element::Zn => 65.38,
            Element::Ga => 69.723,
            Element::Ge => 72.630,
            Element::As => 74.922,
            Element::Se => 78.971,
            Element::Br => 79.904,
            Element::Kr => 83.798,
            Element::Rb => 85.468,
            Element::Sr => 87.62,
            Element::Y => 88.906,
            Element::Zr => 91.224,
            Element::Nb => 92.906,
            Element::Mo => 95.95,
            Element::Tc => 97.,
            Element::Ru => 101.07,
            Element::Rh => 102.91,
            Element::Pd => 106.42,
            Element::Ag => 107.87,
            Element::Cd => 112.41,
            Element::In => 114.82,
            Element::Sn => 118.71,
            Element::Sb => 121.76,
            Element::Te => 127.60,
            Element::I => 126.90,
            Element::Xe => 131.29,
            Element::Cs => 132.91,
            Element::Ba => 137.33,
            Element::La => 138.91,
            Element::Ce => 140.12,
            Element::Pr => 140.91,
            Element::Nd => 144.24,
            Element::Pm => 145.,
            Element::Sm => 150.36,
            Element::Eu => 151.96,
            Element::Gd => 157.25,
            Element::Tb => 158.92535,
            Element::Dy => 162.50,
            Element::Ho => 164.93,
            Element::Er => 167.26,
            Element::Tm => 168.93,
            Element::Yb => 173.05,
            Element::Lu => 174.97,
            Element::Hf => 178.49,
            Element::Ta => 180.95,
            Element::W => 183.84,
            Element::Re => 186.21,
            Element::Os => 190.23,
            Element::Ir => 192.22,
            Element::Pt => 195.08,
            Element::Au => 196.97,
            Element::Hg => 200.59,
            Element::Tl => 204.38,
            Element::Pb => 207.2,
            Element::Bi => 208.98,
            Element::Po => 209.,
            Element::At => 210.,
            Element::Rn => 222.,
            Element::Fr => 223.,
            Element::Ra => 226.,
            Element::Ac => 227.,
            Element::Th => 232.04,
            Element::Pa => 231.04,
            Element::U => 238.03,
            Element::Np => 237.,
            Element::Pu => 244.,
            Element::Am => 243.,
            Element::Cm => 247.,
            Element::Bk => 247.,
            Element::Cf => 251.,
            Element::Es => 252.,
            Element::Fm => 257.,
            Element::Md => 258.,
            Element::No => 259.,
            Element::Lr => 262.,
            Element::Rf => 267.,
            Element::Db => 268.,
            Element::Sg => 269.,
            Element::Bh => 270.,
            Element::Hs => 269.,
            Element::Mt => 277.,
            Element::Ds => 281.,
            Element::Rg => 282.,
            Element::Cn => 285.,
            Element::Nh => 286.,
            Element::Fl => 290.,
            Element::Mc => 290.,
            Element::Lv => 293.,
            Element::Ts => 294.,
            Element::Og => 294.,
        }
    }
}
