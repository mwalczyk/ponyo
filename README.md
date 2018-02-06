# Ponyo
A semi-Lagrangian fluid solver based on Robert Bridson's book, 
_Fluid Simulation for Computer Graphics (2nd Edition)_. This implementation
uses a staggered marker-and-cell (MAC) grid with second-order Runge Kutta 
interpolation during the backwards particle trace. Currently, the pressure
solve is accomplished using the [Gauss-Seidel](https://en.wikipedia.org/wiki/Gauss%E2%80%93Seidel_method) 
method.

<p>
  <img src="https://github.com/mwalczyk/ponyo/blob/master/logo.svg" alt="plume logo" width="100" height="auto"/>
</p>

## Running
Make sure you have [Rust](https://www.rust-lang.org/en-US/install.html) installed. Navigate inside the directory and
execute the command: `cargo build --release`. You can then run the application with `cargo run`. There
are several optional command-line arguments that can be used to configure the simulation:

```
USAGE:
    ponyo.exe [OPTIONS]

FLAGS:
    -h, --help       Prints help information
    -V, --version    Prints version information

OPTIONS:
    -f, --frame_count <COUNT>             Sets the number of frames that will be rendered
    -o, --output_directory <DIRECTORY>    Sets the directory where images will be rendered to
    -r, --resolution <WIDTHxHEIGHT>       Sets the resolution of the simulation, in pixels

```

Note that if you are running the executable via `cargo`, you must place an additional
set of hyphens before any command-line arguments. For example: `cargo run -- -o /renders`.
## Dependencies
Ponyo uses the following crates:
```
image = "*"
pbr = "*"
clap = "*"
```

## References
Throughout the creation of this project, several papers and 
tutorials were particularly helpful:

1. _Fluid Flow for the Rest of Us_ [[link]](https://pdfs.semanticscholar.org/9d47/1060d6c48308abcc98dbed850a39dbfea683.pdf
)
2. _Fluid Simulation for Computer Graphics: A Tutorial in Grid Based
and Particle Based Methods_ [[link]](https://cg.informatik.uni-freiburg.de/intern/seminar/gridFluids_fluid-EulerParticle.pdf
)
3. _Gauss-Seidel Method for Solving Simultaneous Linear Equations_ [[link]](https://www.youtube.com/watch?v=ajJD0Df5CsY)

A few existing repositories were helpful as well, most notably [Incremental
Fluids](https://github.com/tunabrain/incremental-fluids), a series of tutorials 
written by Benedikt Bitterli, whose `projection` method was key to my 
understanding of the pressure solve.