# Sum of Squares Programming for Julia.

| **Documentation** | **Build Status** | **Social** | **References to cite** |
|:-----------------:|:----------------:|:----------:|:----------------------:|
| [![][docs-stable-img]][docs-stable-url] | [![Build Status][build-img]][build-url] | [![Gitter][gitter-img]][gitter-url] | [![DOI][zenodo-img]][zenodo-url] |
| [![][docs-latest-img]][docs-latest-url] | [![Codecov branch][codecov-img]][codecov-url] | [<img src="https://upload.wikimedia.org/wikipedia/en/a/af/Discourse_logo.png" width="64">][discourse-url] | [.bib](https://github.com/jump-dev/SumOfSquares.jl/blob/master/CITATION.bib) |

This packages contains the Sum of Squares reformulation for polynomial optimization.
When used in conjunction with [MultivariatePolynomial](https://github.com/JuliaAlgebra/MultivariatePolynomials.jl) and [PolyJuMP](https://github.com/jump-dev/PolyJuMP.jl), it provides a Sum of Squares Programming extension for [JuMP](https://github.com/jump-dev/JuMP.jl).
Enabling the creation of sum of squares variables and constraints.

## Documentation

- [**STABLE**][docs-stable-url] &mdash; **most recently tagged version of the documentation.**
- [**LATEST**][docs-latest-url] &mdash; *in-development version of the documentation.*

Some presentations on, or using, SumOfSquares (see [here](https://github.com/blegat/SumOfSquaresSlides) for the source code of the presentations):
  * [UPCOMING 13 September 2021 14:00-17:00 CEST] Benoît Legat at [POEMA Learning Week 2](http://poema-network.eu/index.php/news-and-events/project-workshops/13-poema-learning-week-2)
  * [UPCOMING 28 July 2021 16:30-17:00 UTC] Benoît Legat, Marek Kaluba and Tillmann Weisser at [JuMP-dev 2021](https://pretalx.com/juliacon2021/talk/L8DTE3/) [[Slides](https://drive.google.com/file/d/1HtArDFNMQ6IYUqRjSWR3JviJp9xLtSlB/view?usp=sharing)]
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

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-latest-img]: https://img.shields.io/badge/docs-latest-blue.svg
[docs-stable-url]: https://jump.dev/SumOfSquares.jl/stable
[docs-latest-url]: https://jump.dev/SumOfSquares.jl/latest

[build-img]: https://github.com/jump-dev/SumOfSquares.jl/workflows/CI/badge.svg?branch=master
[build-url]: https://github.com/jump-dev/SumOfSquares.jl/actions?query=workflow%3ACI
[codecov-img]: http://codecov.io/github/jump-dev/SumOfSquares.jl/coverage.svg?branch=master
[codecov-url]: http://codecov.io/github/jump-dev/SumOfSquares.jl?branch=master

[gitter-url]: https://gitter.im/JuliaOpt/SumOfSquares.jl?utm_source=share-link&utm_medium=link&utm_campaign=share-link
[gitter-img]: https://badges.gitter.im/JuliaOpt/SumOfSquares.jl.svg
[discourse-url]: https://discourse.julialang.org/c/domain/opt

[zenodo-url]: https://doi.org/10.5281/zenodo.1208672
[zenodo-img]: https://zenodo.org/badge/DOI/10.5281/zenodo.1208672.svg
