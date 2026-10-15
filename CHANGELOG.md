# Changelog

## [0.7.0] - 2024-12-19

### 🚀 Features

- Add script to fix doctests on demand by @kunzaatko


### 🐛 Bug Fixes

- *(docs)* Missing docs by @kunzaatko

- *(CI)* Aqua.jl does not allow differing Project.toml by @kunzaatko


### ⚙️ Miscellaneous Tasks

- *(julia)* [**breaking**] Bump minimum Julia version to 1.8 by @kunzaatko

- Lighten the compat restriction of ColorTypes by @kunzaatko

- *(registrator)* Bump minor version for breaking changes by @kunzaatko


## [0.6.1] - 2024-12-19

### 📚 Documentation

- Improve sampled transfer functions documentation by @kunzaatko


## [0.6.0] - 2024-12-19

### 🚀 Features

- *(core)* [**breaking**] Add apply function for transfer function simulation by @kunzaatko

- *(docs)* Add InterLinks support for external documentation by @kunzaatko


### 🐛 Bug Fixes

- *(docs)* Extra whitespace in documentation causes codeblock by @kunzaatko

- *(docs)* Doc-strings formatting by @kunzaatko

- Test ambiguities with Test by @kunzaatko

- *(tests)* No method for rand with interval by @kunzaatko

- *(CI)* Remove filter from doctest example by @kunzaatko


### 🚜 Refactor

- *(test)* Reorganize test structure and enhance doctests by @kunzaatko


### 📚 Documentation

- Add sampled transfer functions documentation by @kunzaatko


### ⚙️ Miscellaneous Tasks

- Bump minor version by @kunzaatko


## [0.5.0] - 2024-11-18

### 🐛 Bug Fixes

- *(CI)* Bump version of checkout action by @kunzaatko

- LaTeX in bibliography by @kunzaatko

- Add ExactBlackman to documentation by @kunzaatko


### 💼 Other

- Files with dashes instead of underscores by @kunzaatko

- *(tests)* Tests for apodization and edge tapering by @kunzaatko

- *(apodization)* Edgetapering types and documentation by @kunzaatko

- *(tests)* Tests for apodization and edge tapering by @kunzaatko

- *(docs)* Missing documentation and not included docstrings by @kunzaatko

- *(docs)* Rewording documentation by @kunzaatko

- *(tag)* Bump minor version by @kunzaatko


## [0.4.0] - 2024-08-29

### 🚀 Features

- *(MeasuredPSF)* Improve Base.show for MeasuredPSF by @kunzaatko

- *(CI)* Add workflow dispatch as run condition by @kunzaatko

- *(MeasuredPSF)* Infer center as center of the data if not supplied by @kunzaatko

- Accept any dimension by @kunzaatko

- Full refactor of the package architecture by @kunzaatko

- Restructure the package for SampledOTF type by @kunzaatko

- Support non-zero centers for SampledPSF by @kunzaatko

- Add `support`, `overlap` and `cutoff` functions by @kunzaatko


### 🐛 Bug Fixes

- *(MeasuredPSF)* Check that `center` is inbounds by @kunzaatko

- *(CI)* Do not run doctests when building documentation by @kunzaatko

- Tuple(...) -> Tuple for type inference by @kunzaatko

- *(psf)* Normalize the PSF so that sum == 1 by @kunzaatko

- Ambiguous methods for general dimensional psf definitions by @kunzaatko

- Aqua.jl ambiguities by @kunzaatko

- *(docs)* Missing comma in make.jl by @kunzaatko

- *(Project.toml)* Unused dependencies by @kunzaatko

- *(CI)* Compat bounds and stale deps by @kunzaatko

- No longer need Tullio by @kunzaatko

- Missing imported functions by @kunzaatko

- *(CI)* Ambiguities on nightly by @kunzaatko


### 💼 Other

- *(CI)* Add latest release by @kunzaatko

- In show of MOTF/MPSF, dispatch on type instead of using eltype by @kunzaatko

- *(test)* Doctests in runtests.jl by @kunzaatko

- *(tests)* Skip Agua.jl when running locally by @kunzaatko

- *(compat)* Bump compat to v"1.8" due to function `allequal` by @kunzaatko

- *(test)* Doctests only on ubuntu due to unitful show being system dependent by @kunzaatko

- *(tests)* Add tests for Base.interfaces by @kunzaatko

- *(docs)* Update and add refs by @kunzaatko

- *(docs)* Documentation for the new API by @kunzaatko

- Bump compat for Interpolations to 0.15, (keep existing compat)

- `otf_support` function by @kunzaatko

- Method for support and overlap with `img` argument by @kunzaatko

- Argument of otf function by @kunzaatko


### ⚙️ Miscellaneous Tasks

- *(compat)* Update compat entry by @kunzaatko

- *(version)* Tick version for updated API by @kunzaatko


## [0.3.2] - 2023-10-29

### 🚀 Features

- Add possibility to generate for given image by @kunzaatko

- Add possibility to generate `otf_support` for image by @kunzaatko

- Convenience methods to OTF support function by @kunzaatko

- [**breaking**] Return `BitArray` from `otf_support` by @kunzaatko

- Add inclusiveness kwarg to otf_support by @kunzaatko

- Add shift `δ` to `otf` function by @kunzaatko

- Prepare for measured TF implementation by @kunzaatko

- *(CI)* Run CI on PR by @kunzaatko

- *(perf)* Simplify the generation PSF when TF is symmetric by @kunzaatko


### 🐛 Bug Fixes

- Editted the version to be semvered and ignoring research files by @kunzaatko

- [**breaking**] Error types from `ArgumentError` to `DomainError` by @kunzaatko

- Dispatch instead of unions and tests for OTF by @kunzaatko

- Test not written correctly by @kunzaatko

- *(otf)* Error in indexing `Δxy` in the `otf` function by @kunzaatko

- *(CI)* Unsatisfiable requirements for LazyGrids for 1.0 by @kunzaatko

- *(test)* Synchronization doctests by @kunzaatko

- *(test)* Type promotion and doctests by @kunzaatko

- *(docs)* @repl block remove `julia>` by @kunzaatko

- *(docs)* Markdown fixes by @kunzaatko

- ModelTransferFunction not defined before use by @kunzaatko

- *(test)* Use `FourierTools.shift` in testing shift consistency by @kunzaatko

- *(test)* Run tests for Ambiguities for VERSION 1.1 and further by @kunzaatko


### 💼 Other

- Born and Wolf and circular aperture models by @kunzaatko

- *(doc)* `otf` function documentation for otf image generation by @kunzaatko

- Export `ipsf` by @kunzaatko

- *(rename)* ClosedForm{OTF/PSF}ModelTransferFunction -> ClosedForm{OTF/PSF}OTFModel by @kunzaatko

- *(CI)* Coverage processing by @kunzaatko

- Add new compat entry for SpecialFunctions at version 2, (keep existing compat)

- *(docs)* Docs for missing references by @kunzaatko

- `shift` function for SIM by @kunzaatko

- `MeasuredPSF` by @kunzaatko

- Add new compat entry for Unitful at version 1, (keep existing compat)

- *(compat)* Julia ^1.6 by @kunzaatko

- *(Project.toml)* Old syntax for test dependencies by @kunzaatko

- *(rename)* OTFModel -> ModelOTF to confer with MeasuredOTF by @kunzaatko

- *(docs)* Structure and features in docs by @kunzaatko

- Docs structure by @kunzaatko


### 🧪 Testing

- *(Agua.jl)* Aqua.jl automatic code quality tests by @kunzaatko

- *(CI)* Add 1.0 julia version to matrix acording to the compat by @kunzaatko


### ⚙️ Miscellaneous Tasks

- *(rename)* SymmetricPupilFunction -> RadiallySymmetric by @kunzaatko

- *(Agua.jl)* Fix suggestions for Agua.jl by @kunzaatko



