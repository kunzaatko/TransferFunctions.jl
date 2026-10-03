set shell := ["bash", "-cu"]
set dotenv-load

[doc("Default recipe shows available commands")]
[private]
default:
    @just --list

root-dir := justfile_directory()

[doc("Generate and open the package documenation")]
[group("docs")]
[script("fish")]
[working-directory(root-dir)]
docs:
    julia --project=docs --color=yes docs/make.jl
    xdg-open docs/build/index.html

[doc("Run unittests for the package")]
[group("test")]
[script("fish")]
[working-directory(root-dir)]
test:
    julia --project -e "using Pkg; Pkg.test()"
