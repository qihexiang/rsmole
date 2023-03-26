use crate::elements::Element;
use lazy_static::lazy_static;
use map_macro::set;
use ndarray::{arr1, arr2, Array1};
use regex::Regex;
use std::{
    collections::{HashMap, HashSet},
    str::FromStr,
};

pub struct Atom {
    element: Element,
    formal_charge: i8,
    isotope: Option<u8>,
}

impl Atom {
    pub fn get_element(&self) -> Element {
        self.element
    }

    pub fn get_formal_charge(&self) -> i8 {
        self.formal_charge
    }

    pub fn get_isotope_mass(&self) -> f64 {
        self.isotope
            .map_or(self.element.get_average_mass(), |i| i as f64)
    }
}

pub enum ParseError {
    CantAnalyzeLine(String),
    CantAnalyzeToken(String),
    UnknownElement(String),
    InvalidPosition(String),
    DuplicatedUniqId(String),
}

impl FromStr for Atom {
    type Err = ParseError;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        lazy_static! {
            static ref RE: Regex = Regex::new(
                r"(?P<isotope>\d+)?(?P<element>[A-z]([a-z])?)(?P<formal_charge>(+|-)\d+)"
            )
            .unwrap();
        }
        let captured = RE.captures(s);
        if let Some(captured) = captured {
            let element = captured
                .name("element")
                .expect("Always get an element group.")
                .as_str();
            let element = Element::from_str(element)
                .map_err(|_| ParseError::UnknownElement(s.to_string()))?;
            let formal_charge: i8 = if let Some(formal_charge) = captured.name("formal_charge") {
                formal_charge
                    .as_str()
                    .parse()
                    .expect("Matched regex can always parsed as u8")
            } else {
                0
            };
            let isotope: Option<u8> = captured.name("isotope").map(|matched| {
                matched
                    .as_str()
                    .parse()
                    .expect("Always be a unsigned integer")
            });
            Ok(Atom {
                element,
                isotope,
                formal_charge,
            })
        } else {
            Err(ParseError::CantAnalyzeToken(s.to_string()))
        }
    }
}

pub type SelectedGroups = HashMap<String, HashSet<usize>>;

type DescartesMoleItem = (Atom, Array1<f64>, HashSet<String>, String);
pub struct DescartesMole {
    atoms: Vec<DescartesMoleItem>,
}

impl FromStr for DescartesMole {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let lines = if s.contains("\r\n") {
            s.split("\r\n")
        } else {
            s.split("\n")
        }
        .filter(|line| *line != "");
        let mut atoms = vec![];
        for (idx, line) in lines.enumerate() {
            let components = line
                .split(" ")
                .filter(|component| *component != "")
                .collect::<Vec<_>>();
            let atom = components
                .get(0)
                .map(|atom_token| Atom::from_str(*atom_token))
                .unwrap_or(Err(ParseError::CantAnalyzeLine(line.to_string())))?;
            macro_rules! parse_position {
                ($token_idx: expr) => {{
                    components
                        .get($token_idx)
                        .map(|token| {
                            token
                                .parse::<f64>()
                                .map_err(|_| ParseError::InvalidPosition($token_idx.to_string()))
                        })
                        .unwrap_or(Err(ParseError::InvalidPosition(line.to_string())))
                }};
            }
            let x = parse_position!(1)?;
            let y = parse_position!(2)?;
            let z = parse_position!(3)?;
            lazy_static! {
                static ref UID_RE: Regex = Regex::new(r"@([A-Z,a-z])([A-Z,a-z,0-9,_])*").unwrap();
            }
            let (uid, groups) = if components.len() > 4 {
                let (uid, groups) = if UID_RE.is_match(components[4]) {
                    if components.len() > 5 {
                        let groups = &components[5..]
                            .iter()
                            .map(|s| s.to_string())
                            .collect::<HashSet<String>>();
                        (components[4].to_string(), groups.clone())
                    } else {
                        (components[4].to_string(), HashSet::<String>::new())
                    }
                } else {
                    let groups = &components[4..]
                        .iter()
                        .map(|s| s.to_string())
                        .collect::<HashSet<String>>();
                    (format!("@{}", idx + 1), groups.clone())
                };
                (uid, groups)
            } else {
                let groups = HashSet::<String>::new();
                (format!("@{}", idx + 1), groups.clone())
            };
            atoms.push((atom, arr1(&[x, y, z]), groups, uid))
        }
        Ok(DescartesMole { atoms })
    }
}

impl DescartesMole {
    pub fn select(&mut self, group_name: &str, index: usize) -> Option<bool> {
        let (_, _, groups, _) = self.atoms.get_mut(index)?;
        Some(groups.insert(group_name.to_string()))
    }

    pub fn unselect(&mut self, group_name: &str, index: usize) -> Option<bool> {
        let (_, _, groups, _) = self.atoms.get_mut(index)?;
        Some(groups.remove(group_name))
    }

    pub fn get_groups(&self) -> HashSet<String> {
        let mut set = HashSet::new();
        for (_, _, groups, _) in &self.atoms {
            for item in groups.iter() {
                set.insert(item.to_string());
            }
        }
        set
    }

    pub fn get_group_mut(&mut self, group_name: &str) -> Vec<&mut DescartesMoleItem> {
        self.atoms
            .iter_mut()
            .filter(|(_, _, groups, _)| groups.contains(group_name))
            .collect()
    }

    pub fn get_group(&self, group_name: &str) -> Vec<&DescartesMoleItem> {
        self.atoms
            .iter()
            .filter(|(_, _, groups, _)| groups.contains(group_name))
            .collect()
    }

    pub fn remove_atoms(&mut self, group_name: &str) {
        self.atoms
            .retain(|(_, _, groups, _)| !groups.contains(group_name))
    }

    pub fn move_atoms(
        &mut self,
        group_name: &str,
        move_vector: Array1<f64>,
    ) -> Vec<&DescartesMoleItem> {
        let atoms = self.get_group_mut(group_name);
        for (_, position, _, _) in atoms {
            *position = position.clone() + move_vector.clone();
        }
        self.get_group(group_name)
    }

    pub fn rotate_atoms(
        &mut self,
        group_name: &str,
        axis: Array1<f64>,
        radian: f64,
    ) -> Vec<&DescartesMoleItem> {
        let atoms = self.get_group_mut(group_name);
        let (x, y, z) = (axis[0], axis[1], axis[2]);
        let cos = radian.cos();
        let sin = radian.sin();
        let rotate_matrix = arr2(&[
            [
                cos + (1. - cos) * x.powf(2.),
                (1. - cos) * x * y - sin * z,
                (1. - cos) * x * z + sin * y,
            ],
            [
                (1. - cos) * x * y + sin * z,
                cos + (1. - cos) * y.powf(2.),
                (1. - cos) * y * z - sin * x,
            ],
            [
                (1. - cos) * x * z - sin * y,
                (1. - cos) * y * z + sin * x,
                cos + (1. - cos) * z.powf(2.),
            ],
        ]);
        for (_, position, _, _) in atoms {
            *position = position.dot(&rotate_matrix);
        }
        self.get_group(group_name)
    }

    pub fn merge_as_group(
        &mut self,
        group_name: &str,
        structure: DescartesMole,
        keep_origin_groups: bool,
    ) -> Vec<&DescartesMoleItem> {
        if keep_origin_groups {
            for (atom, position, groups, uid) in structure.atoms {
                self.atoms.push((
                    atom,
                    position,
                    {
                        let mut set = set! {};
                        set.extend(groups);
                        set.insert(group_name.to_string());
                        set
                    },
                    uid,
                ))
            }
        } else {
            for (atom, position, _, uid) in structure.atoms {
                self.atoms
                    .push((atom, position, set! {group_name.to_string()}, uid))
            }
        }
        self.get_group(group_name)
    }
}

// fn replace(target: &mut DescartesMole, polymer: DescartesMole, target_prefix: &str, polymer_src: &str, polymer_dst: &str) {
//     let to_replace = target.get_groups();
//     let to_replace = to_replace.iter().filter(|group_name| group_name.starts_with(target_prefix)).collect::<HashSet<_>>();
//     for replace_target in to_replace {
//         let target_src_group_name = replace_target.split(":").collect::<Vec<_>>()[1];
//         let (_, target_src_position, _) = target.get_group(target_src_group_name)[0];
//         let target_dsts = target.get_group(&replace_target);
//         for (_, target_dst_position, _) in target_dsts {
//             let bond_vector = target_src_position - target_dst_position;

//         }
//     }
// }

#[test]
fn overflow_slice() {
    let arr: Vec<usize> = vec![];
    println!("{:?}", &arr[4..])
}
