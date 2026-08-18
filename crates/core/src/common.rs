pub(crate) trait ErrorVecExt<E> {
    fn push_err<T>(&mut self, err: E) -> Option<T>;
}
