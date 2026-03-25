# Contributing to Pavan

## Workflow
1. Fork and clone
2. Create a feature branch
3. Run `make check` (fmt + clippy + test + audit)
4. Add tests and benchmarks for new code
5. Submit PR

## Code Style
- `cargo fmt` enforced
- `cargo clippy -- -D warnings`
- All physics must include correctness tests with known engineering values
- No `unwrap()` or `panic!()` in library code
