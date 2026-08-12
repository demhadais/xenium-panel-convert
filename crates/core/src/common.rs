#[derive(Debug, Clone, Copy)]
#[cfg_attr(feature = "clap", derive(clap::ValueEnum))]
pub enum Species {
    HomoSapiens,
    MusMusculus,
}

#[derive(Debug, Clone, Copy)]
#[cfg_attr(feature = "clap", derive(clap::ValueEnum))]
pub enum Chemistry {
    V1,
    Prime,
}

pub(crate) trait ErrorVecExt<E> {
    fn push_err<T>(&mut self, err: E) -> Option<T>;
}
