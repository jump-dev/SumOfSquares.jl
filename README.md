# SumOfSquares.jl

[![Build Status](https://github.com/jump-dev/SumOfSquares.jl/workflows/CI/badge.svg?branch=master)](https://github.com/jump-dev/SumOfSquares.jl/actions?query=workflow%3ACI)
[![codecov](https://codecov.io/gh/jump-dev/SumOfSquares.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/jump-dev/SumOfSquares.jl)
[![](https://zenodo.org/badge/DOI/10.5281/zenodo.1208672.svg)](https://doi.org/10.5281/zenodo.1208672)

[SumOfSquares.jl](https://github.com/jump-dev/SumOfSquares.jl) is a JuMP
extension that, when used in conjunction with [MultivariatePolynomial](https://github.com/JuliaAlgebra/MultivariatePolynomials.jl)
and [PolyJuMP](https://github.com/jump-dev/PolyJuMP.jl), implements a sum of
squares reformulation for polynomial optimization.

## License

`SumOfSquares.jl` is licensed under the [MIT license](https://github.com/jump-dev/SumOfSquares.jl/blob/master/LICENSE.md).

## Installation

Install `SumOfSquares` using `Pkg.add`:
```julia
import Pkg
Pkg.add("SumOfSquares")
```

## Documentation

See [https://jump.dev/SumOfSquares.jl/stable](https://jump.dev/SumOfSquares.jl/stable)
for the most recently tagged version of the documentation.

See [https://jump.dev/SumOfSquares.jl/dev](https://jump.dev/SumOfSquares.jl/dev)
for the in-development version of the documentation.

## Presentations

Some presentations on, or using, SumOfSquares (see [blegat/SumOfSquaresSlides](https://github.com/blegat/SumOfSquaresSlides)
for the source code of the presentations):

  * Benoît Legat at [JuMP-dev 2023](https://pretalx.com/juliacon2023/talk/XLT8H3/) [[Slides](https://drive.google.com/file/d/1H-_Ot7tP2g7t95r0K_mTaXvGO1ChUPci/view?usp=drive_link)]
  * Benoît Legat, Marek Kaluba and Tillmann Weisser at INFORMS 2022 [[Slides](https://drive.google.com/file/d/1rlsIxgcnKWT436k4MNenjHfgH0UYRLAB/view?usp=share_link)]
  * Benoît Legat at [POEMA Learning Week 2](http://poema-network.eu/index.php/news-and-events/project-workshops/13-poema-learning-week-2)
  * Benoît Legat, Marek Kaluba and Tillmann Weisser at [JuMP-dev 2021](https://pretalx.com/juliacon2021/talk/L8DTE3/) [[Slides](https://drive.google.com/file/d/1HtArDFNMQ6IYUqRjSWR3JviJp9xLtSlB/view?usp=sharing)]
  * Benoît Legat at [INFORMS 2020](http://meetings2.informs.org/wordpress/annual2020/) [[Slides](https://drive.google.com/file/d/1lb8NtOWCikTYm6KRUZCSLYgaUjqIsSyV/view?usp=sharing)]
  * Tillmann Weisser, Benoît Legat, Chris Coey, Lea Kapelevich and Juan Pablo Vielma at [JuliaCon 2019](https://juliacon.org/2019/) [[Slides](https://drive.google.com/open?id=1HiA-praFyejE0Z3nVSpFEv938TAcPjA9)] [[Video](https://www.youtube.com/watch?v=cTmqmPcroFo)]
  * Benoît Legat at [CNLS 2019](https://cnls.lanl.gov/External/showtalksummary.php?selection=7768) [[Slides](https://drive.google.com/open?id=1kNF18C7RY2zi7jcZBMO1PRXtHuvVTFPn)]
  * Benoît Legat at [EURO 2019](https://www.euro2019dublin.com/) [[Slides](https://drive.google.com/open?id=1Wry56NzzL4QBRSwuhP4AlKOe2i2FL7dk)]
  * Benoît Legat at [juliaday Nantes 2019](https://julialang.univ-nantes.fr/programme/) [[Slides](https://drive.google.com/open?id=1pN3G9Pr8jbzK9EEaJ9a6p_qKwSbxb2bo)]
  * Benoît Legat at [Summer School on Numerical Computing in Algebraic Geometry 2018](https://www.mis.mpg.de/calendar/conferences/2018/nc2018.html) [[Poster](https://drive.google.com/open?id=1pf9rdoVEjAnD164rptLki1AG0AH4i88M)]
  * Benoît Legat at [The First Annual JuMP-dev Workshop 2017](https://jump.dev/meetings/mit2017/) [[Slides](https://drive.google.com/file/d/1ea5eSMvMB3jXPuljzNGmMKied-n50YIo/view?usp=sharing)] [[Video](https://youtu.be/kyo72yWYr54)]
  * [Joey Huchette at SIAM Opt 2017](https://docs.google.com/presentation/d/1ASfjB1TdLJmYxT0b6rnyGh9eLbMc-66bTOt3_3yvc90/edit?usp=sharing)

## Citing

See [CITATION.bib](https://github.com/jump-dev/SumOfSquares.jl/blob/master/CITATION.bib).
