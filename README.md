# Ponyo
A semi-Lagrangian fluid solver based on Robert Bridson's book, 
_Fluid Simulation for Computer Graphics (2nd Edition)_.

<p>
  <img src="https://github.com/mwalczyk/ponyo/blob/master/logo.svg" alt="plume logo" width="100" height="auto"/>
</p>

## Running
Make sure you have Rust installed. Navigate inside the directory and
execute the command: `cargo run --release`. This will save frames to 
the directory `images`.

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