# Changelog

## [0.8.0] - 2026-10-15

### 🚀 Features

- *(CI)* Automatic changelog and release on tag push by @kunzaatko

- *(CI)* Use color for the make script in documenter by @kunzaatko

- [**breaking**] Rename IdealCircularPupilOTF to CircularPupilOTF and update docs by @kunzaatko

- Implement bead (microsphere) model by @kunzaatko

- *(extension)* Add Makie extension skeleton by @kunzaatko

- *(ext)* DimensionalData extension and lift compat restrictions by @kunzaatko

- *(CI)* Generate code diagram to the changelog with release of version by @kunzaatko

- Add circulant tensor-based convolution and dependencies by @kunzaatko

- Update CI and Project.toml by @kunzaatko

- Bump version to 0.8.0 by @kunzaatko

- Introduce @require_interface macro and rename `istransferred` to `insupport` by @kunzaatko

- *(docs)* Restructure documentation by @kunzaatko

- Refactor apodization and array types, add interfaces by @kunzaatko

- *(apodization)* Refactor taperedges implementation by @kunzaatko

- *(types)* Add convolution and linear algebra support by @kunzaatko

- *(benchmark)* Introduce benchmarking framework by @kunzaatko

- *(AI)* Add AGENTS.md by @kunzaatko

- *(estimation)* Skeleton for estimation by @kunzaatko

- *(CI)* Test on version names instead of numbers and add `min`, `pre` by @kunzaatko

- *(CI)* Generate figures in the CI by @kunzaatko

- Add Makie scalebar recipe as extension by @kunzaatko

- Make Estimation.bead return a SpatialArray by @kunzaatko

- Functions for positions and `posgrid` for SampledArray by @kunzaatko

- FWHM function for PSF by @kunzaatko

- Better `conv` method for point spread functions by @kunzaatko

- Add `ReflectedArray` type by @kunzaatko

- *(filter)* Add filtering method and some basic tests by @kunzaatko

- *(ext)* Add extension for ImageCore for filtering by @kunzaatko

- Warn when kern does not contain the origin in `kern_padding` by @kunzaatko

- Constants for ReflectedMatrix and ReflectedVector by @kunzaatko

- Add simple contract function circulant-tensors.jl by @kunzaatko

- Add method for inference of padding from Indices by @kunzaatko

- *(extensions)* Add method for `scalebar` to sampled matrices by @kunzaatko

- Introduce comprehensive transfer function models by @kunzaatko

- *(docs)* Restructure documentation and add developer section by @kunzaatko

- *(docs)* Add favicon to documentation by @kunzaatko

- *(plot)* Add default colours for line and text in `scalebar` by @kunzaatko

- Separate the implementation of `posgrid` into `posaxes` and `posgrid` by @kunzaatko

- Use `posaxes` in the plotting of SpatialMatrix by @kunzaatko

- Add `encircled_energy` and `energy_radius` functions by @kunzaatko

- Determine the size of the PSF from the energy radius by @kunzaatko

- Add synthetic data generation module as an extension to Test and Distributions by @kunzaatko

- Add estimation routine for least squares by @kunzaatko

- *(CI)* Remove non-needed steps from the `docs` workflow and update actions by @kunzaatko

- Implement `FFT` module to take care of the `rfft` optimization by @kunzaatko

- *(FFT)* Improve broadcasting system for FFTOut types by @kunzaatko

- *(ext)* Add support for `FFT.fft` on Colorant types by @kunzaatko

- *(test)* Add `setup_makiemaestro!` and `setup_params!` as helper eval functions by @kunzaatko

- Make BorderArray from string/symbol definable with a custom Border type and test by @kunzaatko

- Broadcasting for SampledArrays now return a SampledArray output by @kunzaatko

- Add broadcasting for sampled arrays with check for sampling equality by @kunzaatko

- *(sampling)* Add sample cells and flexible vertex definitions by @kunzaatko

- Add `otf` method by @kunzaatko

- *(docs)* Add column label as a xtick to the documentation of filtering-matrices by @kunzaatko

- Add a `sample_freqs` function on sampled arrays by @kunzaatko

- Add utility functions for printing the models and use them for the printing by @kunzaatko

- Constructor of `OriginAt` with any convertible to `CartesianIndex` and doctests write by @kunzaatko

- Introduce offset-preserving view and axis manipulation utilities by @kunzaatko

- *(psf)* Add Wiener-filter deconvolution by @kunzaatko

- *(release)* Automate changelog and package registration


### 🐛 Bug Fixes

- *(CI)* TagBot workflow update (#12) by @kunzaatko in #12

- *(docs)* Show only used bibliography by @kunzaatko

- *(compat)* [**breaking**] Restrict to julia ^1.9 by @kunzaatko

- *(compat)* Add compat entry for statistics by @kunzaatko

- *(docs)* Add missing doc strings to docs by @kunzaatko

- *(documentation)* Add missing docstrings and references by @kunzaatko

- *(estimation)* Bead function arguments by @kunzaatko

- *(docs)* Reorganize documentation, remove dead links, add missing by @kunzaatko

- *(CI)* Typo by @kunzaatko

- *(docs)* Include missing documents by @kunzaatko

- *(interface)* Add error types and test the correct types and messages for `@require_interface` by @kunzaatko

- *(docs + CI)* Fix the doc tests and included documentation by @kunzaatko

- *(docs)* Fix docstring for Triangular Apodization Function type by @kunzaatko

- *(docs)* Rename doctestssetup to correct file and use in docs/make.jl by @kunzaatko

- *(CI)* Bad `run` field by @kunzaatko

- *(CI)* Remove `lts` version due to `public` keyword use by @kunzaatko

- *(docs)* Correct the external references and missing docstrings by @kunzaatko in #22

- *(CI)* Allow nightly to fail by @kunzaatko

- *(CI)* Separate the nightly tests to allow failure by @kunzaatko

- *(SampledArray)* Base.similar method by @kunzaatko

- *(test)* Do not load MakieMaestro in the doctest setup script by @kunzaatko

- Sampled array constructor type promotion and by @kunzaatko

- Free InterfaceFunctions.jl by @kunzaatko

- *(docs)* Update location of Unitful.jl documentation by @kunzaatko

- *(docs)* Add new types to reference and update the external links by @kunzaatko

- *(CI)* Update the versions of actions by @kunzaatko

- *(CI)* Diagram generation by @kunzaatko

- *(docs)* Warn only for missing docs and fix the errors by @kunzaatko

- *(CI)* Do not run JET tests on nightly in CI by @kunzaatko

- Return nothing from `doctestsetup.jl` to avoid clutter by @kunzaatko

- *(test)* Do not test `JET` on nightly by @kunzaatko

- *(test)* Set the `LINES` and `COLUMNS` for doctests by @kunzaatko

- Rename LinearTransferFunction -> LinearShiftInvariantTransferFunction by @kunzaatko

- *(compat)* Avoid `Fix` for julia 1.11 and earlier and add docstring to `intensity` by @kunzaatko

- *(tests)* Correct the error types thrown and fix the unbound args by @kunzaatko

- Wrong @static condition in version by @kunzaatko

- *(compat)* Version fix for @static and use generic closures in <= v"1.11" by @kunzaatko

- *(docs)* Add missing doc reference to API-reference by @kunzaatko

- *(CI)* Allow the coveralls upload to fail by @kunzaatko

- *(docs)* Add missing references by @kunzaatko

- *(CI)* Remove the unnecessary develop of the package by @kunzaatko

- *(docs)* Markdown format for bulletpoints by @kunzaatko

- Export types in sampled arrays and fix similar Base method by @kunzaatko

- `inferpadding` for indices by @kunzaatko

- *(test)* Update the items that have moved (`bead`, `Beads`) by @kunzaatko

- Allow `conv` with `tf` to be used for colorants by removing the type parameter by @kunzaatko

- *(test)* Fix the failing tests from FFT by @kunzaatko

- *(test)* Documentation of `scalebarformat` by @kunzaatko

- *(test)* Add `BenchmarkTools` to `doctestsetup.jl` by @kunzaatko

- *(test)* Check with JET only on 1.12-rc with function_filter by @kunzaatko

- *(docs)* Include the outputs by @kunzaatko

- *(test)* Doctests setup script by @kunzaatko

- Use CairoMakie in the doctests instead of GLMakie by @kunzaatko

- *(test)* Add the output types of scalebar in doctests by @kunzaatko

- *(ext)* Makie ImageLike requires an Interval returned by @kunzaatko

- *(docs)* Doc tests by @kunzaatko

- *(docs)* Markdown list in GibsonLanni by @kunzaatko

- Add two more methods of BroadcastStyle for disambiguation by @kunzaatko

- Use the eltype for similar in sampled-array by @kunzaatko

- *(test)* Add conditional loading for `doctestsetup` by @kunzaatko

- In `freqgrid` only allow same PixelSize and size by @kunzaatko

- *(docs)* Correct source of `objects.inv` for `DomainSets` by @kunzaatko

- *(build)* Quote success messages in package recipes by @kunzaatko

- *(filter)* Handle mixed real and complex correlation by @kunzaatko


### 💼 Other

- *(CI)* Do not run on 1.12 by @kunzaatko in #20

- *(extension)* MakieMaestro extension by @kunzaatko

- Extension interface by @kunzaatko

- Add `params` method and the constructors from the params for models and augmentations and test by @kunzaatko

- Add just recipes for docs and tests by @kunzaatko

- Add package environment management recipes by @kunzaatko

- Add Julia REPL recipe to justfile by @kunzaatko

- *(deps)* Expand compatible dependency versions by @kunzaatko

- Configure Julia workspace and support Julia 1.13 by @kunzaatko

- *(deps)* Bump julia-actions/cache from 2 to 3 (#29) by @dependabot[bot] in #29

- *(deps)* Bump actions/checkout from 3 to 7 (#28) by @dependabot[bot] in #28


### 🚜 Refactor

- Remove custom NotImplementedError type by @kunzaatko

- File organization by @kunzaatko

- Rename IdealOTFwithCurvature -> CircularPupilOTF by @kunzaatko

- *(estimation)* Bead function and its documentation clean by @kunzaatko

- Significant type restructure by @kunzaatko

- Clean up types, test dependencies and methods by @kunzaatko

- Tests and add methods for util functions by @kunzaatko

- Taperedges and apodize functions by @kunzaatko

- *(arrays)* Refactor circulant and filtering array types by @kunzaatko

- Replace FilteringMatrix indexing functions by properties by @kunzaatko

- *(utils)* Add utility functions + better signatures for existing by @kunzaatko

- *(apodization)* [**breaking**] Rework SineSum and Blackman type parameters and constructors by @kunzaatko

- *(arrays)* Rename OuterInnerArray to Flattened and split type definitions by @kunzaatko

- *(test)* Refactor the tests to include partial runs by @kunzaatko

- *(interfaces)* Use InterfaceFunctions.jl for @interface definitions by @kunzaatko

- Remove ImageCore dependency in favour of extension by @kunzaatko

- Move makie extension to its own folder and add compat for ComputePipeline by @kunzaatko

- *(arrays)* [**breaking**] Refactor the code base and array types by @kunzaatko

- Use alias for SpatialArray{<:Any, 2} in signature by @kunzaatko

- Move the PSF types to a separate directory by @kunzaatko

- Use a util function to determine the `inner_axes` in the tapered array by @kunzaatko

- Separate function for `\sigma` in Gaussian model by @kunzaatko

- *(synthetic-data)* Move synthetic data generation to main module by @kunzaatko

- Allow `psf` to be generated from indices instead of size by @kunzaatko

- Remove unused type parameters in FFT `copyto!` by @kunzaatko

- Add similar for each of the subtype that it is possible to define similar in OffsetArrays with by @kunzaatko

- *(CI)* Rename the `changelog.yml` -> `Changelog.yml` for consistency by @kunzaatko

- Broadcasting of the FFTOut arrays using promotion rules instead of specialized methods by @kunzaatko

- Rename `fftfreqs` to `freqgrid` for consistency by @kunzaatko

- *(borders)* [**breaking**] Rename circular border to periodic by @kunzaatko


### 📚 Documentation

- Fix typo, add Frequency doc, update todo comments, docs style by @kunzaatko

- *(readme)* Add 'made in Czechia' badge by @kunzaatko

- *(comments)* Add the todo comments from files by @kunzaatko

- *(apodization)* Add documentation page with figures for apodization by @kunzaatko

- *(logo)* Add logos to package documentation by @kunzaatko

- Enhance the logos by @kunzaatko

- Add documentation strings and remove unused types by @kunzaatko

- Documentation string for `intensity` of PSF by @kunzaatko

- Add/remove for new/deleted objects documentation by @kunzaatko

- Add documentation for the array types by @kunzaatko

- *(README)* Add logo to README by @kunzaatko

- Specialize the signature for conv and deconv by @kunzaatko

- *(comment)* More explanation on why the doctests are ran in runtests.jl by @kunzaatko

- Add figure to filtering-matrices.md with the column filtering matrix by @kunzaatko

- Add docstrings to conv(!) and corr(!) by @kunzaatko

- Add filter.md doc page by @kunzaatko

- Improve LaTeX rendering and cross-references in docstrings by @kunzaatko

- Add synthetic data doc page by @kunzaatko

- Update API reference and navigation by @kunzaatko

- Update compat of DocumenterCitations by @kunzaatko

- Move the `scalebar!` and `scalebar` docstrings to make them discoverable by @kunzaatko

- *(refactor)* Restructure the filtering docs by @kunzaatko

- Add the missing docs and better phrasing by @kunzaatko

- Round the outputs of the characteristics in airy-disc.md by @kunzaatko

- *(border-arrays)* Document filtering helper functions by @kunzaatko

- *(bibliography)* Normalise BibTeX entries by @kunzaatko

- *(psf)* Remove inaccurate imfilter argument note by @kunzaatko

- *(changelog)* Add historical release notes


### ⚡ Performance

- Remove LazyGrids and stale dependencies in favour of @inline modern julia implementations by @kunzaatko


### 🎨 Styling

- Rename `img` -> `A` and `kern` -> `K` in filtering by @kunzaatko

- Format by @kunzaatko

- Export the types from their respective files where they are defined by @kunzaatko

- Use TF when imported as a nick for TransferFunctions by @kunzaatko

- Format with Julia Language Server by @kunzaatko

- Format with Julia Language Server by @kunzaatko


### 🧪 Testing

- *(compat)* Use CompatHelperLocal for local compat testing by @kunzaatko

- Add tests for interior method by @kunzaatko

- Simplify and reduce the sizes of tests to test faster by @kunzaatko

- Add a test on filtering matrices with ImageCore by @kunzaatko

- Specify the `target_module` in `JET` tests by @kunzaatko

- Add test for correspondence of `filtering_matrix` and `corr` by @kunzaatko

- Add some tests based on coverage misses by @kunzaatko

- Better tests for filtering matrices by @kunzaatko

- Better tests for point spread functions by @kunzaatko

- Better tests for PSF models by @kunzaatko

- Add test for the indices method on the PSF by @kunzaatko

- Add the JET function filter again by @kunzaatko

- Add tests for `setindex!` with `RFFTOut` by @kunzaatko

- Add failing tests for conv in 3D and passing tests for conv in 2D by @kunzaatko

- Update tests for eltype promotion in broadcasting of FFT output types by @kunzaatko

- Add some trivia tests for the filtering by @kunzaatko

- Printing of the models by @kunzaatko


### ⚙️ Miscellaneous Tasks

- *(CompatHelper)* Bump compat for LazyGrids to 1, (keep existing compat) (#13) by @github-actions[bot] in #13

- *(docs)* Add the colprac badge to the repo for conributor instructions by @kunzaatko

- *(git)* Add files to ignore by @kunzaatko

- *(bibliography)* Update bibliography by @kunzaatko

- *(todo)* Remove the unplanned TODO comments by @kunzaatko

- *(git)* Remove globally ignored from gitignore by @kunzaatko

- *(benchmark)* Add scripts for prettifying and version and add git ignore by @kunzaatko

- *(compat)* Update compat and `Project.toml`s by @kunzaatko

- *(CI)* Update the workflow files by @kunzaatko

- *(docs)* Fixes in make.jl and references by @kunzaatko

- *(CI)* Test on nightly by @kunzaatko

- *(compat)* Fix compat for InterfaceFunctions by @kunzaatko

- *(git)* Ignore files that are no longer included in the package by @kunzaatko

- *(rename)* Extensions -> package-extensions by @kunzaatko

- *(format)* Src/types/filtering-matrices.jl by @kunzaatko

- *(docs)* Remove extra blank lines by @kunzaatko

- *(test)* Remove TensorOperations from test dependencies by @kunzaatko

- *(format)* Run Julia format by @kunzaatko

- *(git)* Add cached docs to gitignore by @kunzaatko

- *(compat)* Bump compat of MakieMaestro in docs by @kunzaatko

- *(docs)* Update the bibliography by @kunzaatko

- *(docs)* Update refs by @kunzaatko

- *(CI)* Update the action versions by @kunzaatko

- *(test)* Add download for madril_gray to doctestsetup.jl by @kunzaatko

- *(deps)* Replace CompatHelper with Dependabot by @kunzaatko


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



