---
title: Towards a new Projection Code
author:
  - Martin Ueding (<ueding@hiskp.uni-bonn.de>)
date: 2019-01-12
abstract: >
  The [existing projection
  code](https://github.com/HISKP-LQCD/sLapH-projection) only works for the
  $\rho$ resonance. We want to extend it to multiple particles and arbitrary
  isospin. This note contains the architecture and a survey of possible
  programming languages and libraries for the implementation. Everything that
  is going to be implemented will be written out here.
urlcolor: blue
...

# Specification

The [sLapH contraction code](https://github.com/HISKP-LQCD/sLapH-contractions)
can only work with single particle operators multiplied together. These
operators are just Dirac-$\Gamma$ structure and integer three-momentum $p$. The
spin and isospin projection code needs to take these and give a GEVP for given
total momentum $P^2$ and irrep $\Gamma$.

We have to decide whether and when we want to include baryons. As discussed
with Carsten we will stick with the single cover of the octahedral group and
therefore just mesons for the time being.

The result that we want from the projection is a data frame with these columns:

- Total momentum $P^2$
- Irrep $\Gamma$
- Irrep row $\alpha$
- GEVP row index (source operator ID)
- GEVP column index (sink operator ID)
- HDF5 dataset name
- Complex weighting factor

With this we can take the output of the sLapH contraction code and project into
the irreps.

In order to get this data frame, we will do the first four of these steps:

1. Group theory
2. Isospin projection
3. Spin projection
4. Analytic Wick contractions
5. Actual projection of the correlators from the HDF5 files.

# Architecture

Isospin and spin projection should be independent of each other. The current
projection code does both concepts at the same time and exploits certain
simplifications that come with the $\rho$-channel. In this new code we want to
refrain from any such premature simplifications to use it for any physical
process to come.

## Group theory

We will need to have the following of the octahedral group and each of the
relevant little groups:

- Irreducible representation matrices $D^\Gamma$
- Irreducible continuum matrices $D^J$, the Wigner-D matrices
- Cartesian 3D representation matrices $R_g$

There must be correspondence between the various representations, otherwise we
could not sensibly multiply different representations together.

It still needs to be decided whether we want to have all the little groups
corresponding to every momentum or whether we do rotations to a reference
momentum ourselves.

## Isospin

The isospin projection can be done by hand for the cases we have done so far.
For the $I = 0$, $I_3 = 0$ case for instance we construct a two-pion operator
as
$$
\pi\pi(\vec p_1, \vec p_2)
= \pi^+(\vec p_1) \; \pi^-(\vec p_2)
+ \pi^-(\vec p_1) \; \pi^+(\vec p_2)
- \frac{1}{\sqrt{3}} \pi^0(\vec p_1) \; \pi^0(\vec p_2) \,.
$$

For $I = 1$, $I_3 = 0$ the operator is
$$
\pi\pi(\vec p_1, \vec p_2)
= \frac{1}{\sqrt{2}} \left(\pi^+(\vec p_1) \; \pi^-(\vec p_2)
- \pi^-(\vec p_1) \; \pi^+(\vec p_2) \right) \,.
$$

For $I = 3$, $I_3 = 3$ the operator is just
$$
\pi\pi\pi(\vec p_1, \vec p_2, \vec p_3)
= \pi^+(\vec p_1) \; \pi^+(\vec p_2) \; \pi^+(\vec p_2) \,.
$$

For scalar particles this is very easy. For vector mesons like the $\rho$ we
need to build the various spin states $J_3$.

As pointed out by Marcus this can also be done automatically. As I plan this to
be independent of the remainder, it does not matter so much whether and when we
do this automatically.

## Spin

For the spin part we consult Markus' “projection” chapter and use
Equation (3.34). For the generalization to $N_\text{p}$ particles we
have the following expression:
\begin{align*}
O_\Gamma^{\alpha}(\vec p_\text{cm})^\dagger &=
\sum_{M=-J}^J
\phi_M
\left[
\prod_{i=1}^{N_\text{p}}
\sum_{M_i=-J}^J
\right]
\langle J, M | J_1, M_1, \ldots, J_{N_\text P}, M_{N_\text P} \rangle
\\&\quad\times
\sum_{\beta=1}^{\mathop{\mathrm{nrow}(\Gamma)}}
\tilde\phi_\beta 
\sum_{g \in \mathop{\mathrm{LG}}(\vec p_\text{cm})}
D_{\alpha\beta}^\Gamma(g)^*
\\&\quad\times
\prod_{i=1}^{N_\text{p}}
\sum_{M_i'=-J_i}^{J_i}
\sum_{M_i''=-J_i}^{J_i}
D_{M_i' M_i}^{J_i}(g) \;
D_{M_i'' M_i'}^{J_i}(\tilde g) \;
O_{i M_i''}^{J_i}(R_{\tilde g}^{-1} R_g R_{\tilde g} \vec p_i)^\dagger
\,.
\end{align*}
The variables $J$ and $J_i$ are determined by the physical process that one
investigates.

In order to express this in terms of already executed operator products we need
to rewrite the last line as
\begin{multline*}
\prod_{i=1}^{N_\text{p}}
\sum_{M_i'=-J_i}^{J_i}
\sum_{M_i''=-J_i}^{J_i}
D_{M_i' M_i}^{J_i}(g) \;
D_{M_i'' M_i'}^{J_i}(\tilde g) \;
O_{i M_i''}^{J_i}(R_{\tilde g}^{-1} R_g R_{\tilde g} \vec p_i)^\dagger
\\
\begin{split}
&=
\left(
\prod_{i=1}^{N_\text{p}}
\sum_{\mu_i=-J_i}^{J_i}
\right)
\left[
\prod_{i=1}^{N_\text{p}}
\sum_{M_i'=-J_i}^{J_i}
\sum_{M_i''=-J_i}^{J_i}
D_{M_i' M_i}^{J_i}(g) \;
D_{M_i'' M_i'}^{J_i}(\tilde g) \;
\delta_{M_i'' \mu_i}
\right]
\\&\quad\times
%O_{\mu \nu_\mu}^{J_\mu}(R_{\tilde g}^{-1} R_g R_{\tilde g} \vec p_\mu)^\dagger
O_{\mu_1 \ldots \mu_{N_\text p}}^{J_1\ldots J_{N_\text p}}(R_{\tilde g}^{-1} R_g R_{\tilde g} (\vec p_1, \ldots, \vec p_{N_\text P}))^\dagger
\,.
\end{split}
\end{multline*}

The round parentheses are meant to just limit the scope of the $\prod$-sign
whereas the square brackets are meant to limit the scope of both $\prod$- and
$\sum$-signs. The group element $\tilde g$ is defined as $\vec p_\text{cm} =
R_{\tilde g} \vec p_\text{ref}$. This likely is not unique because we can
always add a rotation around the resulting $\vec p_\text{cm}$ without changing
it. Perhaps this does not matter? Also reflections might need to be excluded.

## Analytic Wick contraction

The analytic Wick contraction can be done with the “Quark Contraction Tool” in
Mathematica. For the $I = 3$, $I_3 = 3$ channel for zero momentum with source
at `x` and sink at `y` we have summands like the following, which we identify
as `C6cD`:

```mathematica
trace[Gamma^5.DE[{up, up}, {x, y}].Gamma^5.DE[{dn, dn}, {y, x}]]^3
```

Also we have the `C6cCD`:

```mathematica
trace[Gamma^5.DE[{up, up}, {x, y}].Gamma^5.DE[{dn, dn}, {y, x}]]
trace[Gamma^5.DE[{dn, dn}, {y, x}].Gamma^5.DE[{up, up}, {x, y}].
      Gamma^5.DE[{dn, dn}, {y, x}].Gamma^5.DE[{up, up}, {x, y}]]
```

And also the `C6cC`:

```mathematica
trace[Gamma^5.DE[{dn, dn}, {y, x}].Gamma^5.DE[{up, up}, {x, y}].
      Gamma^5.DE[{dn, dn}, {y, x}].Gamma^5.DE[{up, up}, {x, y}].
      Gamma^5.DE[{dn, dn}, {y, x}].Gamma^5.DE[{up, up}, {x, y}]]
```

With the sLapH method at zero momentum we have it quite easy as the three
sources and three sinks are indistinguishable. But with general momenta we
cannot simplify it as much. The expression from the contraction code with
sources `x1`, `x2`, `x3` and sinks `x4`, `x5`, `x6` contains `C6cD` diagrams
like this one:

```mathematica
trace[Gamma^5.DE[{up, up}, {x1, x4}].Gamma^5.DE[{dn, dn}, {x4, x1}]]
trace[Gamma^5.DE[{up, up}, {x2, x6}].Gamma^5.DE[{dn, dn}, {x6, x2}]]
trace[Gamma^5.DE[{up, up}, {x3, x5}].Gamma^5.DE[{dn, dn}, {x5, x3}]]
```

But also we have this one here:

```mathematica
trace[Gamma^5.DE[{up, up}, {x1, x4}].Gamma^5.DE[{dn, dn}, {x4, x1}]]
trace[Gamma^5.DE[{up, up}, {x2, x5}].Gamma^5.DE[{dn, dn}, {x5, x2}]]
trace[Gamma^5.DE[{up, up}, {x3, x6}].Gamma^5.DE[{dn, dn}, {x6, x3}]]
```

The difference is just that the exchange $x_5 \leftrightarrow x_6$ was done.

The output from the Quark Contraction Tool is basically a prescription to a
linear combination of diagrams from the contraction code. The positions are
equivalent to momentum labels. From this output we learn the way that we have
to combine the different diagrams to yield a particular isospin.

We have to convert this output into a data frame with the following columns:

- GEVP row index (source operator ID)
- GEVP column index (sink operator ID)
- HDF5 data set name like `C2c_uu_p010.d000.g5_p010.d000.g5`
- Boolean flag whether the correlator needs to be conjugated to compensate for
  reversal of quark flow direction
- Complex weight factor

Ideally we also apply certain symmetry simplifications where we can apply say
the permutation $(13)(24)$ on a `C4cD` diagram and get the same result.

## Spin and isospin

We have prescriptions for both spin and isospin. They both need to be resolved
to yield the wanted linear combination of computed correlation functions.

We have already expressed the spin projection in terms of multi-particle
operators. The isospin projected operators can therefore be resolved there.
From there we can spin project onto given momentum and irrep and so forth. The
result will be multiple-particle operators with definite spin and isospin.
Taking these we can do the Wick contractions.

How do we get the dataset names? Just as we would have to parse the Mathematica
output into HDF5 dataset name *patterns* after doing the Wick contractions
right after the isospin projection, we now need to extract the full dataset
name from the symbolic expression that contains the concrete momenta.

That expression will likely have very many summands. Parsing of it needs to be
done as part of the program. Either with the PyParsing library in Python or
within Mathematica.

# Choice of languages and libraries

We have a natural interface between the group theory part that provides us with
a table of matrices in various irreps and the projection which contains the
physical process.

The projection shall touch the actual numeric correlators as late as possible.
I want to have it compute the spin and isospin projection operator as a data
frame and only then apply it to the data.

For the choice of programming languages and libraries keep in mind that
programming languages are not tools but the very material from which we build
things. They stick around once the project is finished.

## Group theory

For the group theory there are a few options:

Copying tables
  ~ There are various books and papers which list all the things that we need.
  Copying them is just error-prone and we would have to check for consistency
  ourselves. Apparently the `Bethe` Maple library also has it copied from a
  paper.

Maple
  ~ This is bad because neither the university nor the institute has have a
  license for it. We could ask the institute to buy one, though I do not
  believe that the institute can buy the student version for around 150 EUR.
  The `Bethe` library has been used by Markus for the $\rho$-project. Perhaps
  we just ask him to export all we need and hope that we never have to touch it
  again.

Mathematica
  ~ The university has a Mathematica campus license, so although that is
  license riddled as well, at least *we* can use it. The octahedral group is
  not directly available in Mathematica but I managed to get it via the
  representation as a permutation group. I can also construct the little groups
  (also known as stabilizers) but not the irreps.

Sage
  ~ I have been looking into alternatives a bit yesterday and found that Sage
  has the group in some form that one can at least get the multiplication and
  character table, see [this
  snippet](https://sagecell.sagemath.org/?z=eJyNjrFOwzAQhvdIeQdLDLWFY0pYGApLJRgQ4gGqNjrcwzVy4mBfI_z2xIGAQB3QLb_-T9_d2bb3gVgEg8oEf-yjaoGCfW_MGF9sZwldagx2GIBwXxZ5HtRquGU3_9fUOmnnybdW31l0e34tygLGDUM17K7YWXwLxOuxs7nb1WXxpCnzx2nrfT7BN7y6PLfiopbfgfGfpsphK9nmswO5lGwp-dTDVuS_-2A74ov1AQJowsAInh0uxEzyUaVn2kyUi18qJIfptDehk5LvXo8GdGLaQYwY_5ozb774qH8A9juBiQ==&lang=sage).

SciPy
  ~ This has Clebsch-Gordan coefficients and Wigner-3j symbols, but not
  really much more.

We now read in the tables that are generated by Markus.

## Projection operator construction

Also we need Clebsch-Gordan coefficients, but that should be covered in all the
packages that we consider.

There are very likely to be cancellations between different terms. The rotation
matrices contain factors like $1/\sqrt{2}$ and therefore we likely want to do
symbolic algebra. We have these options:

Maple
  ~ Not a good idea as we do not want to depend on it more than strictly
  needed.

Mathematica
  ~ We have a campus license, it is a very powerful language. However it is a
  proprietary tool and our group has limited experience with it. I personally
  do not want to invest too much learning into a proprietary tool.

    Wigner-D matrices are available as `WignerD`.

Sage
  ~ For what we need it will likely be overpowered. As it is a pain to install,
  we should only do that when strictly needed.

SymPy
  ~ This would be my favorite as I and others know Python and it is free
  software.

    Wigner-D matrices are available as `physics.quantum.spin.Rotation`

We chose Mathematica for this as we also do the Wick contractions with
Mathematica. This way we have just one code.

## Wick contraction

Cadabra
  ~ While asking about Wick contractions in SymPy, the author of the
  [Cadabra](https://cadabra.science/) CAS system reached out to me offering
  help implementing Wick contractions in that system. It is programmable in
  Python.

Mathematica
  ~ We have the [Quark Contraction Tool](https://arxiv.org/abs/1603.01576)
  already available.

SymPy
  ~ There is `physics.secondquant.wicks`, which seems to be able to do Wick
  contractions. Also there is `physics.secondquant.contraction` which does a
  contraction between two fermionic operators. I have spend the afternoon of
  2019-01-29 looking into it, but the definition of the fermion creation and
  annihilation operators seem to not support a tensorial structure, at least
  not with the Wick contractions provided. I have asked a [question on Stack
  Overflow](https://stackoverflow.com/q/54426941/653152) about this.

We have chosen Mathematica for this step due to the great “Quark Contraction
Tool”.

## Projecting the data

Once we have the data frame with the desired linear combinations of the
correlators from the contraction code, we just have to load the HDF5 data and
project it. This could be done in various languages.

As HDF5 reading seems to be very slow, we might want to consider adding a
re-packaging step just as Markus does by converting the HDF5 files from the
sLapH contraction code into Pandas HDF5 files.

C++
  ~ We do not really need much for projecting the data, except a HDF5 library.
  Perhaps the performance will be better if we do this in C++. But the most
  likely bottleneck is the HDF5 reading, therefore it does not matter which
  language we use.

Python
  ~ The current projection code uses Python and Pandas. It can work with HDF5
  files.

R
  ~ As the analysis is in R one can consider R for doing the actual projection.
  We have HDF5 routines there as well. The advantage would be that the people
  who do not know Python could work on this part of the code.

# Design and implementation

Most of the code will be done in Mathematica. This section contains the design
and implementation of the different tasks along the way.

As the Wolfram Language is very functional and I know some Haskell, I made
great use of this. The following table gives you pointers to some of the
operators used and their correspondence in other languages.

| Long | Infix | Python | R | C++ |
| --- | --- | --- | --- | --- |
| [`f[x]`](http://reference.wolfram.com/language/ref/Prefix.html) | `f @ x` or `x // f` | `f(x)` | `f(x)` | `f(x)` |
| [`Map[f, x]`](http://reference.wolfram.com/language/ref/Map.html) | `f /@ x` | `map(f, x)` | `lapply(x, f)` | `std::transform` |
| [`Composition[f, g]`](http://reference.wolfram.com/language/ref/Composition.html) | `f @* g` | — | — | — |
| [`Function[x, x^2]`](http://reference.wolfram.com/language/ref/Function.html) | `#^2 &` | `lambda` | `function` | `[](…){…}` |
| [`Part[xs, i]`](http://reference.wolfram.com/language/ref/Part.html) | `xs[[i]]` | `xs[i-1]` | `xs[[i]]` | `xs[i-1]` |
| [`Association[a -> b]`](http://reference.wolfram.com/language/ref/Association.html) | `<| a -> b |>` | `{a: b}` | `list(a = b)` | [`std::map`](https://en.cppreference.com/w/cpp/container/map) |
| [`Dot[a, b]`](http://reference.wolfram.com/language/ref/Dot.html) | `a . b` | `a @ b` | `a %*% b` | — |

Common patterns in functional programming are [*tacit
programming*](https://en.wikipedia.org/wiki/Tacit_programming) and expressing
loops in terms of list comprehensions or *maps*. I use some of the former and
plenty of the latter. This might make the code harder to read for a programmer
who is not familiar with the Wolfram Language.

The functions from my `sLapHProjection` package will be explained along the way
and typeset in bold monospace font such that they stand out and can be found
easier (but there's always Ctrl-F).

## Group theory

### Data format for group elements and irreps

In `sciebo/Lattice/Representations` we have a bunch of text files that contain
the group elements and all the irreducible representations for the full group
and subgroup.

In the file `Oh-elements.txt` we find a CSV table with elements like the following:

| Operator |  $\alpha/\pi$ | $\beta/\pi$ | $\gamma/\pi$ |
| :--- | ---: | ---: | ---: |
| `E` | 0 | 0 | 0 |
| `C2x` | 0 | 1 | 1 |
| `C2y` | 0 | 1 | 0 |
| `C2z` | 0 | 0 | 1 |
| `C31+` | 0 | 0.5 | 0.5 |

These are the group elements and their action parametrized in terms of the
three Euler angles $\alpha$, $\beta$ and $\gamma$.

Then there is `Little-group-elements.txt`, a TSV table that contains the subset
of operators for each subgroup associated with a momentum vector $d$. The
entries look like this:

| $\vec{d}$ | $\mathrm{LG}(\vec{d})$ | Operator | $\alpha/\pi$ | $\beta/\pi$ | $\gamma/\pi$ |
| --- | --- | --- | ---: | ---: | ---: |
| (0, 0, 1) | `C4v` | `E` | 0 | 0 | 0 |
| (0, 0, 1) | `C4v` | `C4z+` | 0 | 0 | $0.5$ |
| (0, 0, 1) | `C4v` | `C4z-` | 0 | 0 | $-0.5$ | | (0, 0, 1) | C4v | C2z | 0 | 0 | $1$ |
| (0, 0, 1) | `C4v` | `sigma_x` | 0 | 1 | $1$ |
| (0, 0, 1) | `C4v` | `sigma_y` | 0 | 1 | $0$ |
| (0, 0, 1) | `C4v` | `sigma_d1` | 0 | 1 | $0.5$ |
| (0, 0, 1) | `C4v` | `sigma_d2` | 0 | 1 | $-0.5$ |

The representations for the individual subgroups are in files like
`C2v-(-1,-1,0)-representations.txt`. There also is
`Oh-(0,0,0)-representations.txt` which contains all group elements because the
little group is the group itself:

| Irrep | Oh | Lg | $\alpha/\pi$ | $\beta/\pi$ | $\gamma/\pi$ | row | col | matrix element |
| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `A1g` | `E` | `E` | 0 | 0 | 0 | 1 | 1 | 1. |
| `A1g` | `C2x` | `C2x` | 0 | 1 | 1 | 1 | 1 | 1. |
| `A1g` | `C2y` | `C2y` | 0 | 1 | 0 | 1 | 1 | 1. |
| `A1g` | `C2z` | `C2z` | 0 | 0 | 1 | 1 | 1 | 1. |
| `Eg` | `E` | `E` | 0 | 0 | 0 | 1 | 1 | 1. |
| `Eg` | `E` | `E` | 0 | 0 | 0 | 1 | 2 | 0. |
| `Eg` | `E` | `E` | 0 | 0 | 0 | 2 | 1 | 0. |
| `Eg` | `E` | `E` | 0 | 0 | 0 | 2 | 2 | 1. |
| `Eg` | `C2x` | `C2x` | 0 | 1 | 1 | 1 | 1 | 1. |
| `Eg` | `C2x` | `C2x` | 0 | 1 | 1 | 1 | 2 | 0. |
| `Eg` | `C2x` | `C2x` | 0 | 1 | 1 | 2 | 1 | 0. |
| `Eg` | `C2x` | `C2x` | 0 | 1 | 1 | 2 | 2 | 1. |

The “Irrep” column denotes the irrep name, “Oh” denotes the name of the group
element as part of the $O_h$ group and “Lg” the group element name it has
within the little group. Little groups belonging to the momenta that are
related by a lattice rotation are equivalent to each other but contain
different members of the original group.

### Reading the lattice irreps

In my Mathematica package there is the function **`ReadDataframe`** which reads
such a text file. The function **`ConvertMatrixElement`** will convert that
matrix element column into Mathematica expressions. Use **`ReadIrreps`** for a
convenient composition.

The irrep matrices are represented in the *long format*. For instance that
`C2x` group element in the `Eg` irrep is a $2 \times 2$ matrix that is written
out as
$$
D^{E^+}(C_{2x}) = \begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix} \,.
$$

There are two ways that one can represent tensors:

1.  A high-dimensional array of the actual values. Each dimension corresponds
    to one index. Contractions are performed by looping over all indices that
    stay and then sum up the elements along the contracted dimension.

    The tracking of the index labels it a bit tedious. This also does not scale
    well when indices are added later as the dimensionality of the arrays has to
    be increased. Without proper labelling, this will quickly spiral out of
    control.

    Memory usage is good, though.

2.  A data frame in the *long data format*. Contractions are done by *group-by*
    and *summarize* operations. To contract an index one groups by all the
    indices that are to stay and summarizes the sum of the values.

    This scheme easily scales with arbitrary many indices and adding indices
    later on is not a problem. Memory usage is higher, though zeros can be
    omitted easily.

    It seems that Mathematica has the notion of
    [`Dataset`](http://reference.wolfram.com/language/ref/Dataset.html), which
    means that one can just `Import` a CSV or TSV file and has a data frame as
    with R or Pandas. Also the `Dataset` supports having a complex matrix as
    column type, which would allow for a hybrid approach if we wanted to.

3.  There is another way, using nested
    [`Association`](http://reference.wolfram.com/language/ref/Association.html)
    instances. They are just like the `dict` in Python and the `std::map` in
    C++. This way we could make associations that map from some named group
    element to an association that maps from indices to a c-number.

    As associations can be called like functions, we can really nicely work
    with them. So let's make up a nonsense association like this:

    ```mathematica
    assoc = Association[foo -> bar, 1 -> "one", x^2 -> Association]
    ```

    When you now call `assoc[x^2]`, you get the function `Association`. Or you
    call it with `assoc[1]` and get `"one"`. Missing keys are some sort of
    error condition via say `Missing["KeyAbsent", bar]`.

    Markus was immediately sold on this approach, and I quite like it as well.

I have implemented Item 3. The $D_{\alpha\beta}^\Gamma(g)$ are represented as
an association of this form:
$$ \vec d \to \Gamma \to g \to (\alpha, \beta) \to D_{\alpha\beta}^\Gamma(g) \,. $$
You can think of this as a function taking a momentum vector $\vec d$ and
giving you a new function which accepts an irrep $\Gamma$ and gives you a
function ….

### Creating the Cartesian representation

In the file `Oh-elements.txt` we just have the names of the group elements and
the Euler angles. We read those in with **`ReadEulerAngles`** and obtain a
handy association of the form
$$ g \to (\alpha, \beta, \gamma) \,.$$

The Cartesian representation comes for free, there is
[`EulerMatrix`](http://reference.wolfram.com/language/ref/EulerMatrix.html)
which just takes such a vector. This will give us the mapping $g \to R_g$ that
we need. Great, we're done here if the definition of the Euler angles matches.

### Creating the spin representation

The spin-$J$ representations are for free as well, we just use
[`WignerD`](http://reference.wolfram.com/language/ref/WignerD.html) which is a
mapping
$$ (j, m_1, m_2, \psi, \theta, \phi) \to D^J(g) \,. $$

### Clebsch-Gordan coefficients

The Clebsch-Gordan coefficients are implemented as
[`ClebschGordan`](http://reference.wolfram.com/language/ref/ClebschGordan.html),
so this is shootin' fish in a barrel as well.

# Documentation of functions in the Mathematica module

The following is a list of the functions defined in the Mathematica module
`sLapHProjection`. From what I gathered the documentation mechanism in
Mathematica is not as formal as say Doxygen in C++ or Rd in R. This and the
need for a bunch of math expressions made me chose an external documentation.

These functions are ordered by name to make finding stuff easy at the expense
of a coherent reading from top to bottom.

-   $\mathtt{EulerGTilde}(\vec p_\text{cm})$

    Finds a set of Euler angles $(\psi, \theta, \phi)$ such that the rotation
    matrix coverts the reference momentum into center-of-mass momentum, i.e.
    $$\vec p_\text{cm} = R(\psi, \theta, \phi) \; \vec p_\text{ref} \,.$$ This
    transformation is not unique, but that *should* not be a problem as we sum
    over the little group elements anyway.

-   $\mathtt{MomentumRef}(\vec p_\text{cm})$

    Gives the reference momentum $\vec p_\text{ref}$ for given center-of-mass
    momentum $\vec p_\text{cm}$.
