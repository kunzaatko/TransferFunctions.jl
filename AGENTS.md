## AGENTS.md

### Test:

The test command is:
`julia --project -e 'using Pkg; Pkg.test()'`

When adding a new feature, it is best to create a separate file for the feature's tests within the `test/` directory.
This helps to keep the `runtests.jl` file organized and maintainable.

Tests should comply with the structure of `test/runtests.jl`. Add new features to test independently, similar to the
existing `"obligatory"` and `"optional"` sections.

Example:
```julia
@testset "MyNewFeature" begin
    # Test code here
end
# or if it is part of a larger test suite for some part of the package
@cond_testset "MyNewFeature" begin
    include("my_new_feature.jl")
end
```

### Documentation:

Documentation is written in Markdown using the [`Documenter.jl`](https://github.com/JuliaDocs/Documenter.jl) markdown
syntax. To build the documentation, use the following command from the project root:
`cd docs/ && julia --project --color=yes make.jl && cd ..`

Documentation should be brief and focused. Advanced features can be covered in reference sections. Examples should be
included directly in markdown example blocks (e.g., ````@example feature ... ````), not via `include()`. Examples should
use real-life scenarios and `Base` Julia types, concentrating on the features of this package. Object names in examples
should be descriptive (e.g., `FlyingVehicle`, `fuel_consumed`).

Documentation structure preferences:
- The `index.md` file should include a short, brief summary of the package's features as a list.
- The `quickstart.md` section should provide a concise introduction to enable the user to quickly start using the
  package main features.
- The `manual.md` section should contain more detailed examples, explanations and should touch on the philosophy and the
  reasoning of the package / API structure within the text of the various feature explanations and demonstrations. The
  math and/or logic behind any calculations and algorithms should be included where necessary. It should however remain
  to the point and simple enough for a novice user to understand.

Documentation structure preferences:
- The `index.md` file should include a short, brief summary of the package's features.
- The `quickstart.md` section should provide a concise introduction.
- The `manual.md` section should contain more detailed examples and explanations.
