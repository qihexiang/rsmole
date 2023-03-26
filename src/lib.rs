use std::ops::Add;

pub mod elements;
pub mod molecule;

pub fn add<T>(left: T, right: T) -> T
where T: Add<Output = T>
{
    left + right
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn it_works() {
        let result = add(2, 2);
        assert_eq!(result, 4);
    }
}
