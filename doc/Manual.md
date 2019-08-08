---
 # vim: spell tw=79

 # Document meta data:
title: New Projection Code
author:
  - Martin Ueding (<ueding@hiskp.uni-bonn.de>)
date: 2019-01-12
abstract: >
  The [existing projection
  code](https://github.com/HISKP-LQCD/sLapH-projection) only works for the
  $\rho$ resonance. We want to extend it to multiple particles and arbitrary
  isospin. This note contains the architecture and a survey of possible
  programming languages and libraries for the implementation. Also the design
  and implementation is covered here, including documentation for all functions
  in the Mathematica package.

 # LaTeX presentation options:
urlcolor: blue
toc: true
numbersections: true
documentclass: scrreprt

 #mainfont: Libertinus Serif
sansfont: Noto Sans
 #monofont: Noto Mono
 #mathfont: Libertinus Math
...

# Specification

The [sLapH contraction code](https://github.com/HISKP-LQCD/sLapH-contractions)
can only work with single particle operators multiplied together. These
operators are just Dirac-$\Gamma$ structure and integer three-momentum $p$. The
spin and isospin projection code needs to take these and give a correlator
matrix for given total momentum $P^2$ and irrep $\Gamma$.

We have to decide whether and when we want to include baryons. As discussed
with Carsten we will stick with the single cover of the octahedral group and
therefore just mesons for the time being.

The result that we want from the projection is a data frame with these columns:

- Total momentum $P^2$
- Irrep $\Gamma$
- Irrep row $\alpha$
- Correlator matrix row index (source operator ID)
- Correlator matrix column index (sink operator ID)
- HDF5 dataset name
- Boolean indicating whether the diagram needs to be complex conjugated
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
\sum_{M_i=-{J_i}}^{J_i}
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

- Correlator matrix row index (source operator ID)
- Correlator matrix column index (sink operator ID)
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

GAP

  ~ A commenter on my [Stack Exchange
    question](https://mathematica.stackexchange.com/q/188767/1507) mentioned
    [GAP](https://www.gap-system.org/). There is seems to be rather easy to
    generate the groups. It is built in or can be created from the generators.
    Then it has the character tables and the conjugacy classes available right
    away.

    ```
    Oh := SmallGroup(48,48);
    Oh := Group( (1,3,4,2)(5,7,8,6), (2,5,3)(4,6,7), (1,8)(2,7)(3,6)(4,5) );
    Display(CharacterTable(Oh));
    Display(ConjugacyClasses(Oh));
    ```

Maple

  ~ This is bad because neither the university nor the institute has have a
    license for it. We could ask the institute to buy one, though I do not
    believe that the institute can buy the student version for around 150 EUR.
    The `Bethe` library has been used by Markus for the $\rho$-project. Perhaps
    we just ask him to export all we need and hope that we never have to touch
    it again.

Mathematica

  ~ The university has a Mathematica campus license, so although that is
    license riddled as well, at least *we* can use it.
  
    The octahedral group is not directly available in Mathematica but I managed
    to get it via the representation as a permutation group. I can also
    construct the little groups (also known as stabilizers) but not the irreps.

    There also is the
    [GTPack](https://www.frontiersin.org/articles/10.3389/fphy.2018.00086/full)
    library which might add what we need. We have not looked into this in
    detail yet.

Sage

  ~ I have been looking into alternatives a bit yesterday and found that Sage
    has the group in some form that one can at least get the multiplication and
    character table, see [this
    snippet](https://sagecell.sagemath.org/?z=eJyNjrFOwzAQhvdIeQdLDLWFY0pYGApLJRgQ4gGqNjrcwzVy4mBfI_z2xIGAQB3QLb_-T9_d2bb3gVgEg8oEf-yjaoGCfW_MGF9sZwldagx2GIBwXxZ5HtRquGU3_9fUOmnnybdW31l0e34tygLGDUM17K7YWXwLxOuxs7nb1WXxpCnzx2nrfT7BN7y6PLfiopbfgfGfpsphK9nmswO5lGwp-dTDVuS_-2A74ov1AQJowsAInh0uxEzyUaVn2kyUi18qJIfptDehk5LvXo8GdGLaQYwY_5ozb774qH8A9juBiQ==&lang=sage):

    ```python
    import sage.groups.matrix_gps.finitely_generated

    K.<v> = sage.groups.matrix_gps.finitely_generated.CyclotomicField(8)
    a = v - v^3
    i = v^2
    Octa = MatrixGroup([(-1+i)/2, (-1+i)/2, (1+i)/2, (-1-i)/2],
                       [(1+i)/a, 0, 0, (1-i)/a])

    print('Character table')
    print(Octa.character_table())
    ```

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
    Python. I have described this project to the author and he told me that he
    was going to try it a bit, but I never heard back from him.

Mathematica

  ~ We have the [Quark Contraction Tool](https://arxiv.org/abs/1603.01576)
    already available.

SymPy

  ~ There is `physics.secondquant.wicks`, which seems to be able to do Wick
    contractions. Also there is `physics.secondquant.contraction` which does a
    contraction between two fermionic operators.

    I have spend the afternoon of 2019-01-29 looking into it, but the
    definition of the fermion creation and annihilation operators seem to not
    support a tensorial structure, at least not with the Wick contractions
    provided. I have asked a [question on Stack
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
    language we use. Declaring the HDF5 data types in C or C++ is rather
    bothersome. We can of course just copy them from the contraction code, but
    having them dynamically is likely easier.

Python

  ~ The current projection code uses Python and Pandas. It can work with HDF5
    files. I do like Python as a language, but knowledge of that language is
    not widespread in our work group.

R

  ~ As the analysis is in R one can consider R for doing the actual projection.
    We have HDF5 routines there as well. The advantage would be that the people
    who do not know Python could work on this part of the code.

I choose R for the implementation. This will also allow bundling it closer with
the Lüscher-Analysis package.

# Design and implementation (Mathematica)

Most of the code will be done in Mathematica. This section contains the design
and implementation of the different tasks along the way.

As the Wolfram Language is very functional and I know some Haskell, I made
great use of this. The following table gives you pointers to some of the
operators used and their correspondence in other languages. Mathematica
functions are in general hyperlinked directly to the official documentation.
You can also just type `?x` to get the help for `x`.

| Long | Infix | Python | R |
| --- | --- | --- | --- | --- |
| `f[x]` | `f @ x` or `x // f` | `f(x)` | `f(x)` |
| [`Apply`]`[f, a]` | `f @@ a` | `f(*a)`[^kwargs] | `do.call(f, a)` |
| [`Map`]`[f, x]` | `f /@ x` | `map(f, x)` | `lapply(x, f)` |
| [`Composition`]`[f, g]` | `f @* g` | —[^pycomp] | [`purrr::compose`](https://purrr.tidyverse.org/reference/compose.html) |
| [`Function`]`[x, x^2]` | `#^2 &` | `lambda` | `function` |
| [`Part`]`[xs, i]` | `xs[[i]]` | `xs[i-1]` | `xs[[i]]` |
| [`Rule`] | `->` | — | — |
| [`RuleDelayed`] | `:>` | — | — |
| [`Association`]`[a -> b]` | `<| a -> b |>` | `{a: b}` | `list(a = b)`[^rlist] |
| [`Dot`]`[a, b]` | `a . b` | `a @ b` | `a %*% b` |
| [`NonCommutativeMultiply`] | `**` | — | — |
| [`ReplaceAll`] | `/.` | — | — |
| [`StringJoin`] | `<>` | `+` | `paste0` |

[^kwargs]:

    In the Wolfram Language the [`List`] instance `a` could also contain
    [`Rule`] instances, therefore making a function call with positional and
    named arguments possible. R's `do.call` also uses `names(a)` as the names
    for named arguments. In Python one cannot mix a list and a dictionary,
    therefore one would have to call `f(*a, **b)` where `b` was a dictionary
    with string keys containing the named arguments.

[^pycomp]:

    Not directly available, but [there are
    ways](https://stackoverflow.com/a/24047214/653152).

[^rlist]:

    This is not entirely correct, though. In R, the names of the list have to
    be strings. For a given `list` instance `l`, we find that
    `typeof(names(l))` is `character` (homogenious string vector). This means
    that we can only have a mapping from strings to arbitrary types, but not
    from arbitrary types to arbitrary types. The Wolfram Language does not have
    this limitation. Also Python only requires the dictionary key to be
    hashable, otherwise there are no constraints.

Note that in R one can create user-defined infix operators, so writing the
following would allow writing R code similar to infix Wolfram Language code:

```r
`%.%` <- purrr::compose
`%/@%` <- function (f, l) lapply(l, f)
`%@@%` <- do.call
```

Common patterns in functional programming are [*tacit
programming*](https://en.wikipedia.org/wiki/Tacit_programming) and expressing
loops in terms of list comprehensions or *maps*. I use some of the former and
plenty of the latter. This might make the code harder to read for a programmer
who is not familiar with the Wolfram Language.

From what I gathered the built-in documentation mechanism in Mathematica is not
as formal as say Doxygen in C++ or Rd in R. The formal one which is used for
the official documentation is too much overhead for my little package.
Therefore I chose to explain the functions from my `sLapHProjection` package
towards the end of each section in this document.

The notation of the function arguments will be mixed text and math. Sometimes
it is more expressive to give a mathematical expression instead of just the
name that the argument has in code. Optional parameters can be recognized by
the default value behind a colon, which is the Wolfram Language notation for
default values.

Since the Wolfram Language is a functional one, the most straightforward way to
define a constant is as a function with zero arguments. In Haskell there is no
distinction between constants and argument less functions at all.

You need to properly install the Quark Contraction Tool. For this download the
`qct.m` from their [GitHub repository](https://github.com/djukanovic/qct) and
then follow [Wolfram's installation instructions](http://support.wolfram.com/kb/5648).

---

-   **`MonitoredMap`**(function, list, label)

    A version of [`Map`] that just maps the *function* over the *list*. It
    shows a progress bar together with the given *label*. With sensible labels
    it can therefore be nested.

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

The irrep matrices are represented in the *long format*. For instance that
`C2x` group element in the `Eg` irrep is a $2 \times 2$ matrix that is written
out as
$$ D^{E^+}(C_{2x}) = \begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix} \,. $$

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

    It seems that Mathematica has the notion of [`Dataset`], which means that
    one can just `Import` a CSV or TSV file and has a data frame as with R or
    Pandas. Also the `Dataset` supports having a complex matrix as column type,
    which would allow for a hybrid approach if we wanted to.

3.  There is another way, using nested [`Association`] instances. They are just
    heterogeneous mappings from one thing to another. This corresponds to the
    `dict` in Python. The `std::map` in C++ is similar, but requires the same
    type among the keys and the same time among the values. Using this for our
    implementation, we could make associations that map from some named group
    element to an association that maps from indices to a c-number.

    As associations can be called like functions, we can really nicely work
    with them. So let's make up a nonsense association like this:

    ```mathematica
    assoc = Association[foo -> bar, 1 -> "one", x^2 -> Association]
    ```

    When you now call `assoc[x^2]`, you get the function `Association`. Or you
    call it with `assoc[1]` and get `"one"`. Missing keys like `bar` are some
    sort of error condition via say `Missing["KeyAbsent", bar]`.

Markus was immediately sold on this approach 3, and I quite like it as well.
This is now implemented.

---

-   **`ReadDataframe`**(filename)

    Reads a CSV formatted table with a header and generates a `Dataset` from
    it.

-   **`ConvertMatrixElement`**(dataset)

    The matrix elements of the irreducible representations,
    $D^\Gamma_{\alpha\beta}(g)$ are stored as strings in the CSV file. Using
    [`ToExpression`] these are parsed into Mathematica symbols. This function
    takes the dataset and replaces the “matrix element” column with the parsed
    one.

-   **`ReadIrreps`**(filename)

    The composition of `ConvertMatrixElement` and `ReadDataframe`.

-   **`MatrixElementsToAssociation`**(dataset)

    Converts the “row”, “col” and “matrix element” into mappings of the form
    $(\alpha, \beta) \to D^\Gamma_{\alpha\beta}(g)$. The result is a
    [`Dataset`] with just one row per irrep and little group element and an
    [`Association`] which contains all the matrix elements.

-   **`IrrepsToAssociation`**(matrixElementsAssociations)

    Continuing with the result from `MatrixElementsToAssociation` this function
    then builds a more deeply nested [`Association`] which has the form
    $$ \Gamma \to g \to (\alpha, \beta) \to D_{\alpha\beta}^\Gamma(g) \,. $$

-   **`DatasetToAssociations`**(dataset)

    Composition of `IrrepsToAssociation` and `MatrixElementsToAssociation`.

-   **`ExtractMomentumFromFilename`**(filename)

    Given a filename like `C2v-(0,1,1)-representations.txt` it extracts the
    three-momentum vector $\vec p_\text{cm}$.

-   **`IrrepDGammaAssoc`**()

    The $D_{\alpha\beta}^\Gamma(g)$ are represented as an association of this
    form:
    $$ \vec d \to \Gamma \to g \to (\alpha, \beta) \to D_{\alpha\beta}^\Gamma(g) \,. $$
    You can think of this as a function taking a momentum vector $\vec d$ and
    giving you a new function which accepts an irrep $\Gamma$ and gives you a
    function which ….

### Creating the Cartesian representation

In the file `Oh-elements.txt` we just have the names of the group elements and
the Euler angles.

The Cartesian representation comes for free, there is [`EulerMatrix`] which
just takes such a vector. This will give us the mapping $g \to R_g$ that we
need. Great, we're done here if the definition of the Euler angles matches.

---

-   **`ReadEulerAngles`**(filename)

    Reads the file with Euler angles and produces a handy association of the
    form
    $$ g \to (\alpha, \beta, \gamma) \,.$$

-   **`EulerAnglesAssoc`**()

    All the Euler angles already read in available as a global constant.

-   **`MomentumRefScalar`**($|\vec p_\text{cm}|^2$)

    Gives the reference momentum. These are hand-chosen to be:
    $$
    \begin{aligned}
    0 &\to (0, 0, 0) \\
    1 &\to (0, 0, 1) \\
    2 &\to (1, 1, 0) \\
    3 &\to (1, 1, 1) \\
    4 &\to (0, 0, 2) \,.
    \end{aligned}
    $$

    This only works only up to $|\vec p_\text{cm}|^2 \leq 8$ because after that
    the mapping is not unique any more, one can write 9 as $9 = 0^2 + 0^2 +
    3^2$ but also as $9 = 2^2 + 2^2 + 1^2$.

-   **`MomentumRef`**($\vec p_\text{cm}$)

    Gives the reference momentum $\vec p_\text{ref}$ for given center-of-mass
    momentum $\vec p_\text{cm}$. As this function is implemented via
    `MomentumRefScalar`, the same restrictions apply.

-   **`EulerGTilde`**($\vec p_\text{cm}$)

    Finds a set of Euler angles $(\psi, \theta, \phi)$ such that the rotation
    matrix coverts the reference momentum into center-of-mass momentum, i.e.
    $$\vec p_\text{cm} = R(\psi, \theta, \phi) \; \vec p_\text{ref} \,.$$ This
    transformation is not unique, but that *should* not be a problem as we sum
    over the little group elements anyway.

-   **`MatrixRGTilde`**($\vec p_\text{cm}$)

    Simply computes the rotation matrix corresponding to the Euler angles from
    \texttt{EulerGTilde}.

-   **`MomentumTransform`**($\vec p$, $\vec p_\text{cm}$, $(\psi, \theta, \phi)$)

    The momentum transformation that is needed in the argument of the
    single-particle operator:
    $$ R_{\tilde g}^{-1} R_g R_{\tilde g} \,. $$

### Creating the spin representation

The spin-$J$ representations are for free as well, we just use [`WignerD`]
which is a mapping
$$ (j, m_1, m_2, \psi, \theta, \phi) \to D^J(g) \,. $$

### Clebsch-Gordan coefficients

The Clebsch-Gordan coefficients are implemented as [`ClebschGordan`]. We do
need to implement higher Clebsch-Gordan coefficients ourselves, though.

This is the relation that I came up with by figuring that we first couple two
spins and then couple that to the next one.
\begin{multline*}
\langle J, M | j_1, m_1, j_2, m_2, j_3, m_3, \ldots \rangle
\\=
\sum_{\tilde J = |j_2 - j_3|}^{j_2 + j_3}
\sum_{\tilde M = - \tilde J}^{\tilde J}
\langle J, M | j_1, m_1, \tilde J, \tilde M \rangle
\langle \tilde J, \tilde M | j_2, m_2, j_3, m_3, \ldots \rangle \,.
\end{multline*}
If the “$\ldots$” are empty, then we have the usual Clebsch-Gordan coefficients
that are already available. This is the stop of the recursion.

Creating a recursion relation is not easy at is seems. I have [asked on Physics
Stack
Exchange](https://physics.stackexchange.com/questions/459408/clebsch-gordan-coefficients-for-more-than-2-particles),
and in one answer I was told that this is correct though not unique.

I was pointed to the [Racah $W$
coefficient](https://en.wikipedia.org/wiki/Racah_W-coefficient). These do not
seem to be implemented in Mathematica, though.

For coupling pions these factors will all be 1 anyway, so we could consider not
bothering with them for just now.

---

-   **`HigherClebschGordan`**($\{ j_i \}$, $\{ m_i \}$, $(J, M)$)

    Computes the Clebsch-Gordan coefficient for coupling all the spins $(j_i,
    m_i)$ together to $(J, M)$. Note that the API is different from
    [`ClebschGordan`] as here the $j_i$ and $m_i$ are grouped among each other
    instead of grouped by $i$.

    In contrast to [`ClebschGordan`] this function does not emit any warnings
    in case the coupling is unphysical.

## Spin

We implement spin independent of isospin here. The spin projection will just
create an expression that contains `SingleOperator`, which can then later be
replaced with the actual operators used. But isospin projection gives us one
big multi-particle operator. I let it compute with the single operators and
then massage the expression such that the products of operators are directly
adjacent. Then I use a pattern replacement to replace them all with one
multi-particle operator.

Markus has chosen the phase vector $\tilde \phi_\beta$ as $\delta_{\beta 1}$,
effectively fixing the column of the irrep. In my implementation I just expose
it as a parameter to the user and remove the sum over the $\beta$.

The formula is implemented from the right to left. We have `MakeSingleOperator`
applying the Wigner-$D$-matrices to the single operators. Then
`MakeMultiOperator` creates a product of the single operators. `MakeGroupSum`
sums over these operator products, `MakeMagneticSum` sums over the $M$ and
$M_i$ quantum numbers. The last step is not finished as the higher
Clebsch-Gordan coefficients are not finished.

From the spin projected operator we extract the actual momenta

---

-   **`MakeSingleOperator`**($\vec p_i$, $\vec p_\text{cm}$, $\vec \Psi_g$,
    $J_i$, $M_i$, $i$)

    Creates an appropriate linear combination of the single particle operator
    $O$ from
    $$
    \sum_{M_i'=-J_i}^{J_i}
    \sum_{M_i''=-J_i}^{J_i}
    D_{M_i' M_i}^{J_i}(g) \;
    D_{M_i'' M_i'}^{J_i}(\tilde g) \;
    O_{i M_i''}^{J_i}(R_{\tilde g}^{-1} R_g R_{\tilde g} \vec p_i)^\dagger \,.
    $$
    The group element $g$ is specified via its Euler angles $\vec \Psi_g$.

    For instance taking $\vec p_i = (1, 1, 1)$, $\vec p_\text{cm} = (1, 1, 1)$,
    $\vec \Psi_g = (\pi/2, 1, 1)$, $J_i = 1$, $M_i = 1$, $i = 1$ gives us the
    following expression:

    ```mathematica
    I ConjugateTranspose[SingleOperator[1, 1, -1, {1, 1, -1}]]
    ```

    The function `SingleOperator` is just a placeholder defined next.

    We can also take a more interesting case where the Wigner $D^J$ matrices
    are not so simple:

    ```mathematica
    MakeSingleOperator[{0, 0, 1}, {1, 1, 1}, Pi {1/2, 1/2, 0}, 2, 0, 1]
    ```

    We then get the following linear combination of operators:

    ```mathematica
    -(1/2) Sqrt[3/2]
       ConjugateTranspose[SingleOperator[1, 2, -2, {0, 1, 0}]] - 
     1/2 ConjugateTranspose[SingleOperator[1, 2, 0, {0, 1, 0}]] - 
     1/2 Sqrt[3/2] ConjugateTranspose[SingleOperator[1, 2, 2, {0, 1, 0}]]
   ```

-   **`SingleOperator`**($i$, $J_i$, $M_i$, $\vec p$)

    Just a placeholder without a definition. It holds the single particle
    operator $O_{i M_i}^{J_i}(\vec p)$. These will later be replaced with the
    isospin projected operators.

-   **`MakeMultiOperator`**($\{ \vec p_i \}$, $\vec \Psi_g$, $\{ J_i \}$, $\{
    M_i \}$)

    Creates the multi particle operator as a product of single particle
    operators. The center-of-mass momentum $\vec p_\text{cm}$ is taken to be
    the sum of the individual momenta $\vec p_i$.

    It corresponds to this expression:

    $$
    \prod_{i=1}^{N_\text{p}}
    \mathtt{SingleOperator}(i, J_i, M_i, \vec p) \,.
    $$

    For example we use the following code:

    ```mathematica
    MakeMultiOperator[{{1, 1, 1} - {1, 0, 0}, {1, 0, 0}}, {1, 0, Pi},
        {1, 1}, {1, 0}]
    ```

    This then gives us the following product of two single operators:

    ```
    (-E^I ConjugateTranspose[
        SingleOperator[1, 1, 1, {Sin[1], -Cos[1], 1}]]) ** 
     ConjugateTranspose[SingleOperator[2, 1, 0, {-Cos[1], -Sin[1], 0}]]
    ```

-   **`MakeGroupSum`**($\Gamma$, $\alpha$, $\beta$, $\{ \vec p_i \}$, $\{ J_i \}$, $\{ M_i \}$)

    This computes the sum over the group elements:
    $$
    \sum_{g \in \mathop{\mathrm{LG}}(\vec p_\text{cm})}
    D_{\alpha\beta}^\Gamma(g)^* \;
    \mathtt{MakeMultiOperator}(\{ \vec p_i \}, \vec \Psi_g, \{ J_i \}, \{ M_i \})
    $$

    The center-of-mass momentum is computed as the sum of the individual momenta.

-   **`MakeMagneticSum`**($\Gamma$, $\alpha$, $\beta$, $\{ \vec p_i \}$, $\{ J_i \}$, $\phi$)

    This function performs the sum over all the magnetic quantum numbers $M$
    and $\{ M_i \}$. The resulting expression is
    \begin{align*}
    O_\Gamma^{J\alpha\beta}(\vec p_\text{cm})^\dagger &=
    \sum_{M=-J}^J
    \phi_M
    \left[
    \prod_{i=1}^{N_\text{p}}
    \sum_{M_i=-J_i}^{J_i}
    \right]
    \langle J, M | J_1, M_1, \ldots, J_{N_\text P}, M_{N_\text P} \rangle
    \\&\quad\times
    \texttt{MakeGroupSum}(\Gamma, \alpha, \beta, \{ \vec p_i \}, \{ J_i \}, \{ M_i \}) \,.
    \end{align*}

-   **`MomentumToString`**($\vec p$)

    Converts an integer momentum vector like `{1, 0, -1}` into the string
    `"10-1"`.

-   **`ExtractMomenta`**(expr)

    Extracts all the momenta from the `SingleOperator`s used. The resulting
    expression contains summand of this form:

    ```mathematica```
    (-(1/2) + (I Sqrt[3])/2) DTMomenta[{1, 0, 1}, {0, 1, 0}]
    ```

    By replacing the `DTMomenta` with something else, one can extract just the
    momenta from this. One has to be a bit careful as Mathematica will try to
    distribute the common factor onto the elements in the lists if one replaces
    it with something simple like [`List`]. Therefore we use something which
    Mathematica does not know and therefore does not try to simplify.

-   **`momentaToRules`**(momenta, location) --- [not exported]

    Takes a product of momentum [`List`]s in a [`NonCommutativeMultiply`] and
    converts them into an [`Association`]. *location* shall be either `"so"` or
    `"si"`

    Calling this function with


    ```mathematica
    MomentaToRules[{0, 1, 1} ** {1, 0, 0}, "si"]
    ```

    gives:

    ```mathematica
    <|"psi1" -> "011", "psi2" -> "100"|>
    ```

-   **`DTMomenta`**

    An undefined symbol which acts as a data type. Used with [`ReplaceAll`].

-   **`MomentaToAssoc`**($\{ \vec p_i \}$, location)

    Converts the momentum terms with [`Association`]s containing terms like

    ```mathematica
    DTMomenta[{0, 1, 1}, {1, 0, 0}]
    ```

    to an association like this:

    ```mathematica
    DTMomentaAssoc[<|"psi1" -> "011", "psi2" -> "100"|>]
    ```

    The name is build up from `p`, then the *location* and consecutive numbers.
    It is wrapped in `DTMomentaAssoc` such that Mathematica does not try to
    distribute numeric factors to the arguments.

-   **`DTMomentaAssoc`**

    An undefined symbol which acts as a data type. Used with [`ReplaceAll`].

-   **`MomentaToAssocSourceSink`**(source, sink)

    Takes two expressions containing `DTMomenta` for the source and sink. They
    get converted into `DTMomentaAssoc` expressions like the following:

    ```mathematica
    DTMomentaAssoc[<|"psi1" -> "011", "psi2" -> "100",
        "pso1" -> "011", "pso2" -> "100"|>]
    ```

## Isospin

Isospin and the Wick contractions are independent of the spin. Therefore we can
do the contractions on the isospin operators and then pattern match that to the
spin expressions. This allows us to go a bit further before connecting isospin
and spin.

We define a $\pi^+$ operator that has two spin, one color and one
position/momentum index as follows:

```mathematica
\[Pi]Plus[s1_, s2_, c1_, x1_] :=
    -FieldB[up, c1, s1, x1] ** (Gamma^5)[SI[{s1, s2}]] **
    Field[dn, c1, s2, x1];
```

Similarly we define the other two pion operators.

From these we build for instance an $I = 2$ operator.

```mathematica
\[Pi]\[Pi]I2[s1_, s2_, s3_, s4_, c1_, c2_, x1_, x2_] :=
    \[Pi]Plus[s1, s2, c1, x1] ** \[Pi]Plus[s3, s4, c2, x2]; 
```

For reasons not entirely clear to me we need to define the “bar” version of
that operator ourselves. This likely is because the QCT does not assume
anything about the Dirac algebra behind the `Gamma` and therefore cannot really
do anything about that. Anyway, the second operator is this:

```mathematica
\[Pi]\[Pi]I2Bar[s1_, s2_, s3_, s4_, c1_, c2_, x1_, x2_] :=
    \[Pi]Minus[s3, s4, c2, x2] ** \[Pi]Minus[s1, s2, c1, x1]; 
```

A few pion operators are defined as part of the library such that they do not
have to be created for every projection projection from scratch.

---

-   **`\[Pi]Plus`**($s_1$, $s_2$, $c_1$, $x_1$)

    $$ [\pi^+]_{s_1 s_2}^{c_1}(x_1) =
    - \bar u_{s_1}^{c_1}(x_1) \; \gamma^5_{s_1 s_2} \; d_{s_2}^{c_1}(x_1) $$

-   **`\[Pi]Minus`**($s_1$, $s_2$, $c_1$, $x_1$)

    $$ [\pi^+]_{s_1 s_2}^{c_1}(x_1) =
    \bar d_{s_1}^{c_1}(x_1) \; \gamma^5_{s_1 s_2} \; u_{s_2}^{c_1}(x_1) $$

-   **`\[Pi]Zero`**($s_1$, $s_2$, $c_1$, $x_1$)

    $$ [\pi^0]_{s_1 s_2}^{c_1}(x_1) = \frac{1}{\sqrt 2} \left(
    \bar u_{s_1}^{c_1}(x_1) \; \gamma^5_{s_1 s_2} \; u_{s_2}^{c_1}(x_1)
    - \bar d_{s_1}^{c_1}(x_1) \; \gamma^5_{s_1 s_2} \; d_{s_2}^{c_1}(x_1)
    \right) $$

-   **`\[Pi]\[Pi]I1`**($s_1$, $s_2$, $s_3$, $s_4$, $c_1$, $c_2$, $x_1$, $x_2$)

    $$ [\pi\pi^{I=1}]_{s_1 s_2 s_3 s_4}^{c_1 c_2}(x_1, x_2) =
    [\pi^+]_{s_1 s_2}^{c_1}(x_1) \; [\pi^-]_{s_3 s_4}^{c_2}(x_2)
    - [\pi^+]_{s_3 s_4}^{c_2}(x_2) \; [\pi^-]_{s_1 s_2}^{c_1}(x_1) $$

-   **`\[Pi]\[Pi]I1Bar`**($s_1$, $s_2$, $s_3$, $s_4$, $c_1$, $c_2$, $x_1$, $x_2$)

    $$ [\bar{\pi\pi}^{I=1}]_{s_1 s_2 s_3 s_4}^{c_1 c_2}(x_1, x_2) =
    [\pi^-]_{s_3 s_4}^{c_2}(x_2) \; [\pi^+]_{s_1 s_2}^{c_1}(x_1)
    - [\pi^-]_{s_1 s_2}^{c_1}(x_1) \; [\pi^+]_{s_3 s_4}^{c_2}(x_2) $$

-   **`\[Pi]\[Pi]I2`**($s_1$, $s_2$, $s_3$, $s_4$, $c_1$, $c_2$, $x_1$, $x_2$)

    $$ [\pi\pi^{I=2}]_{s_1 s_2 s_3 s_4}^{c_1 c_2}(x_1, x_2) =
    [\pi^+]_{s_1 s_2}^{c_1}(x_1) \; [\pi^+]_{s_3 s_4}^{c_2}(x_2) $$

-   **`\[Pi]\[Pi]I2Bar`**($s_1$, $s_2$, $s_3$, $s_4$, $c_1$, $c_2$, $x_1$, $x_2$)

    $$ [\bar{\pi\pi}^{I=2}]_{s_1 s_2 s_3 s_4}^{c_1 c_2}(x_1, x_2) =
    [\pi^-]_{s_3 s_4}^{c_2}(x_2) \; [\pi^-]_{s_1 s_2}^{c_1}(x_1) $$

## Wick contractions

With the appropriate isospin operators defined, we can perform the Wick
contraction using the “quark contraction tool”:

```mathematica
wc = WickContract[\[Pi]\[Pi]I2Bar[s1, s2, s3, s4, c1, c2, x1, x2] **
    \[Pi]\[Pi]I2[s5, s6, s7, s8, c5, c6, x3, x4]];
```

We can then perform the quark contractions on this:

```mathematica
qc = QuarkContract @ wc
```

When we take a look at the output (using `TraditionalForm`) we get something
that looks like this:[^vimregex]
\begin{align*}
&\text{tr}\left(\gamma^5\;S^{\text{up}}(p_4,p_1)\;\gamma^5\;S^{\text{dn}}( p_1,p_4)\right) \text{tr}\left(\gamma^5\;S^{\text{up}}(p_3,p_2)\;\gamma^5\;S^{\text{dn }}(p_2,p_3)\right)
\\&\quad
+\text{tr}\left(\gamma^5\;S^{\text{up}}(p_3,p_1)\;\gamma^5\;S^{\text{dn}}(p_1,p_3)\right) \text{tr}\left(\gamma^5\;S^{\text{up}}(p_4,p_2)\;\gamma^5\;S^{\text{dn }}(p_2,p_4)\right)
\\&\quad
-\text{tr}\left(\gamma^5\;S^{\text{dn}}(p_2,p_3)\;\gamma^5\;S^{\text{up}}(p_4,p_2)\;\gamma^5\;S^{\text{dn}}(p_1, p_4)\;\gamma^5\;S^{\text{up}}(p_3,p_1)\right)
\\&\quad
-\text{tr}\left(\gamma^5\;S^{\text{up}}(p_4,p_1)\;\gamma^5\;S^{\text{dn}}(p_2,p_4)\;\gamma^5\;S^{\text{up}}(p_3,p_2)\;\gamma^5\;S^{\text{dn}}(p_1,p_3)\right)
\end{align*}

[^vimregex]:

    Unfortunately Mathematica does not give this output directly like that, it
    needs `TeXForm` and these replacements in Vim to get it as above:

    ```vim
    s/\\text{Gamma}/\\gamma/g
    s/trace/tr/g
    s/\\text{x\(\d\)}/p_\1/g
    s/\./\\;/g
    ```

The result is something that we will recognize as “C4cD” and “C4cC”. If we are
at zero momentum and there is no distinction between them, this can be
simplified down further to
$$ 2 \cdot \text{C4cD} - 2 \cdot \text{C4cC} \,. $$
In this general case we cannot do that, though.

### Normalizing trace expressions

In the contraction code in
[`src/DiagramSpec.cpp`](https://github.com/HISKP-LQCD/sLapH-contractions/blob/8433f69825f5d5b097f6f30c045af037b06997b9/src/DiagramSpec.cpp)
we have documentation describing the diagrams that can be computed:

Box:
$$
C_\text{C4cB} =
\langle
\Gamma_\mathtt{Op0}
D_\mathtt{Q1}^{-1}(t|t') \Gamma_\mathtt{Op1}
\gamma_5 D_\mathtt{Q2}^{-1}(t'|t')^\dagger \gamma_5
\Gamma_\mathtt{Op2}
D_\mathtt{Q3}^{-1}(t'|t)
\Gamma_\mathtt{Op3}
\gamma_5 D_\mathtt{Q0}^{-1}(t|t)^\dagger \gamma_5
\rangle
$$
Cross:
$$
C_\text{C4cC} =
\langle
\Gamma_\mathtt{Op0}
D_\mathtt{Q1}^{-1}(t|t')
\Gamma_\mathtt{Op1}
\gamma_5 D_\mathtt{Q2}^{-1}(t|t')^\dagger \gamma_5
\Gamma_\mathtt{Op2}
D_\mathtt{Q3}^{-1}(t|t')
\Gamma_\mathtt{Op3}
\gamma_5 D_\mathtt{Q0}^{-1}(t|t')^\dagger \gamma_5
\rangle
$$
Direct:
$$
C_\text{C4cD} =
\langle
\Gamma_\mathtt{Op0}
D_\mathtt{Q1}^{-1}(t|t')
\Gamma_\mathtt{Op1}
\gamma_5 D_\mathtt{Q0}^{-1}(t|t')^\dagger \gamma_5
\rangle
\langle
\Gamma_\mathtt{Op2}
D_\mathtt{Q3}^{-1}(t|t')
\Gamma_\mathtt{Op3}
\gamma_5 D_\mathtt{Q2}^{-1}(t|t')^\dagger \gamma_5
\rangle
$$
Vacuum:
$$
C_\text{C4cV} =
\langle
\Gamma_\mathtt{Op0}
D_\mathtt{Q1}^{-1}(t|t)
\Gamma_\mathtt{Op1}
\gamma_5 D_\mathtt{Q0}^{-1}(t|t)^\dagger \gamma_5
\rangle
\langle
\Gamma_\mathtt{Op2}
D_\mathtt{Q3}^{-1}(t'|t')
\Gamma_\mathtt{Op3}
\gamma_5 D_\mathtt{Q2}^{-1}(t'|t')^\dagger \gamma_5
\rangle
$$

In order to translate the expressions that we get from the Quark Contraction
Tool we need to massage them until they match a diagram from the contraction
code. There are three things that need to be done:

-   The trace is cyclic. Therefore similar expressions can be equivalent. The
    trace expression used, `trace` from `qct`, does not let Mathematica know
    that it is cyclic. Therefore we manually cycle them until the up-type
    propagator starting from the lowest vertex (to be defined) is the second
    element. The first element shall be the $\Gamma$-structure. The first
    vertex is the source vertex with the lowest number. If there is no source
    vertex, then the lowest sink vertex is to be used.

    The key observation is that [expressions can be used like
    lists](https://reference.wolfram.com/language/tutorial/PartsOfExpressions.html).
    With this we can access the different parts of expressions. The arguments
    of a function are subexpressions.

-   The direction of quark flow might be in the wrong direction. By adding a
    transposition to the argument of the trace one can reverse all the
    propagators. This might change the value of the trace if the Dirac
    structures give some signs under transposition.

-   In all the charged/conjugated diagrams the position of the down-type
    propagators is fixed. If from the Wick contraction we happen to have them
    in other positions we need to exchange all flavors in order to get to the
    form that is implemented. This is done by applying a conjugate transpose to
    the argument of the trace. The resulting diagram will need to be
    conjugated.

---

-   **`RotateGammaToFront`**(expression)

    Takes an expression like 

    ```mathematica
    Gamma^5.DE[{"dn", "dn"}, {so[2], si[1]}].
        Gamma^5.DE[{"up", "up"}, {si[2], so[2]}]
    ```

    and cyclicly permutes the expression if it does not start with a `Gamma`.
    The resulting expression always starts with `Gamma`.

-   **`Starts`**(traceContent)

    Gives the starting times of the propagators in the trace content
    expression. The result could look like this:

    ```mathematica
    {si[1], so[2], si[2], so[1]}
    ```

-   **`IndexOfFirst`**(traceContent)

    Determines which propagator should be the first one in the expression. If
    the result is $n$, then the whole expression should be rotated left by $2(n
    - 1)$ iterations.

-   **`NormalizeTrace`**(traceContent)

    Cyclicly permutes the trace such that the Dirac structure corresponding to
    the first vertex is at the front.

-   **`NormalizeTraceRecursive`**(expression)

    Recursively traverses an expression and applies `NormalizeTrace` on every
    `trace` expression encountered.

    The recursion can encounter three cases:

    - The length of the expression is zero, it is just a simple number or so.
      In that case the value is just the expression given.

    - The expresssion matches `trace[_]`, where `_` is a placeholder. In this
      case the argument is run through `NormalizeTrace`.

    - Else this function is mapped over all parts of the expression.

-   **`ReplacePropagators`**(expression)

    Replaces the QCT propagators with a simpler expression. For some reason the
    propagators contain the flavor twice redundantly and also contain sink and
    source time. For the later replacement with diagrams it is sufficient to
    have the flavor once and the source time only. This function does this
    replacement:

    ```mathematica
    qct`DE[{f_, f_}, {_, t_}] :> prop[f, t]
    ```

-   **`FlowReversalRules`**()

    A list of rules that perform reversal of quark flow by taking the transpose
    of the content of the trace. We would need to know what the transpose of
    the Dirac structure is, which depends on the basis, which is not
    implemented yet. Therefore this is limited to the case of $\gamma_5$ which
    will be sufficient as long as we scatter just pions.

-   **`FlavorSwitchRules`**()

    A list of rules that exchange the flavors and might also change the
    direction of quark flow. Similar to the `FlowReversalRules` this currently
    only supports the $\gamma_5$ Dirac structure.

### Contractions into dataset name templates

At the stage of the Wick contractions we still know which labels correspond to
source and sink vertices. We can track this information by adding the undefined
functions `so` and `si` which label the source and sink vertices instead of
`x1` and so on. The contractions are then defined with operators evaluated like
so:

```mathematica
wc = WickContract[
  \[Pi]\[Pi]I2Bar[s1, s2, s3, s4, c1, c2, so[1], so[2]] **
  \[Pi]\[Pi]I2[s5, s6, s7, s8, c5, c6, si[1], si[2]]];
```

In the contracted form we have propagator expressions like these:

```mathematica
prop["up", so[2]]
```

These have to be mapped to HDF5 dataset names. Actually they have to be mapped
to templates ([`StringTemplate`]) that can then be used with [`TemplateApply`]
to insert the momenta. The following replacement can be used to replace an
expression with a `C4cD` diagram:


```mathematica
qct`trace[qct`Gamma^g1_ . prop["up", so[so1_]].
  qct`Gamma^g2_ . prop["dn", si[si2_]]]
qct`trace[qct`Gamma^g3_ . prop["up", so[so3_]].
  qct`Gamma^g4_ . prop["dn", si[si4_]]] :> 
TemplateApply[
  "C4cD_uuuu_" <> MakeTemplate[4],
  <|"g1" -> g1, "g2" -> g2, "g3" -> g3, "g4" -> g4,
    "x1" -> "`pso" <> ToString @ so1 <> "`",
    "x2" -> "`psi" <> ToString @ si2 <> "`",
    "x3" -> "`pso" <> ToString @ so3 <> "`",
    "x4" -> "`psi" <> ToString @ si4 <> "`"|>],
```

The resulting expression after applying these rules is a sum of strings, which
is exactly what we want. In a later stage we can insert the momenta
corresponding to the spin projected operators.

---

-   **`MakeTemplate`**($n$)

    Creates a template for an $n$-particle correlation function consisting of
    expressions like

    ```mathematica
    "g`g1`.p`x1`.d000"
    ```

    separated by underscores and with increasing numbers for `g1` and `p1`.

-   **`DatasetNameRules`**()

    List of delayed rules ([`RuleDelayed`]) that can be used with with
    [`ReplaceAll`] to replace propagator expression with HDF5 dataset name
    templates ([`StringTemplate`]). The replacements are expressions like
    these:
    
    ```mathematica
    "C2cD_uu_g5.p`pso1`.d000_g5.p`psi1`.d000"
    ```

    The following mappings are performed:


## Wick contraction and spin

We now have the needed linear combination of contraction code diagrams. And
using `MomentaToAssocSourceSink` we have a prescription for inserting the
momenta. Inserting every set of momenta into every term of the diagram
expression gives us the desired result, a linear combination of actual HDF5
dataset names.

---

-   **`CombineIsospinAndSpin`**(correlator templates, momenta associations)

    Takes an expression of *correlator templates* containing strings like

    ```mathematica
    "C4cD_uuuu_g5.p`pso1`.d000_g5.p`psi1`.d000_g5.p`pso2`.d000_g5.p`psi2`.d000"
    ```

    and an expression containing *momenta associations* like 

    ```mathematica
    DTMomentaAssoc[<|"pso1" -> "011", "pso2" -> "100", "psi1" -> "011",
        "psi2" -> "100"|>]
    ```

    and gives an expression containing strings like these:

    ```mathematica
    "C4cD_uuuu_g5.p001.d000_g5.p011.d000_g5.p110.d000_g5.p100.d000"
    ```

    This expression then contains everything that is needed to project actual
    data.

### Export to data frame

We want to export the expression in a form like this:

| datasetname | conjugate | re | im |
| --- | --- | ---:| ---: |
| `C4cB_uuuu_g5.p001.d000_…` | False | 1. | 0. |
| `C4cB_uuuu_g5.p001.d000_…` | False | -0.5 | -0.8660254037844386 |
| `C4cB_uuuu_g5.p001.d000_…` | False | 0.5 | -0.8660254037844386 |
| `C4cB_uuuu_g5.p001.d000_…` | False | -0.5 | 0.8660254037844386 |

In CSV format it will be this:

```csv
datasetname,conjugate,re,im
"C4cB_uuuu_g5.p001.d000_…","False",1.,0.
"C4cB_uuuu_g5.p001.d000_…","False",-0.5,-0.8660254037844386
"C4cB_uuuu_g5.p001.d000_…","False",0.5,-0.8660254037844386
"C4cB_uuuu_g5.p001.d000_…","False",-0.5,0.8660254037844386
```

---

-   **`StringExpressionToAssociation`**(expression)

    Takes an expression containing a linear combination of strings and converts
    them into an association where the strings are the keys and the c-number
    valued factors are the values.

    It uses [`Coefficient`] to extract the leading coefficient in front of the
    strings. Also it prefixes the string `conj:` in front of diagram names that
    need to be conjugated.

    This function is based on a [a post by
    Kuba](https://mathematica.stackexchange.com/a/191718/1507) with an implicit
    MIT license.

-   **`NeedsConjugation`**(name)

    Checks whether the *name* starts with `conj:` and returns the list (bare
    name, True/False), where “bare name” is the name without the `conj:`.

-   **`DatasetnameAssocToCSV`**(association, filename)

    Takes the *association* which maps HDF5 dataset names to complex numbers
    (like the result from `StringExpressionToAssociation`) and converts that
    into a CSV structure with these columns:

    - `datasetname`
    - `conjugate`
    - `re`
    - `im`

    This function is based on a [a post by
    Kuba](https://mathematica.stackexchange.com/a/191718/1507) with an implicit
    MIT license.

## Creating a full correlator matrix

### Interface to numerical code

So far we only have functions that generate a list of HDF5 dataset names that
will yield a single correlation function. We of course want to build a full
correlator matrix (with multiple correlators) and then many of these correlator
matrices. This means that we will have the following independent variables
describing the correlators:

- All total momenta $\vec d$ for given $\vec d^2$ (and not just total momentum
  $\vec d^2$)
- Irrep
- Irrep rows
- Correlator matrix row and column (operator ids & relative momenta $\vec q_i$)

The nature of the operator and the relative momenta shall be passed down into
the analysis. This can be done with the meta data notation that Markus has
implemented in his projection code. Although it would be nice to use the same
parser in R, writing it out in a machine readable format would be even better.
Therefore I favor writing out the meta data in the YAML format.

The information per correlator (`datasetname`, `conjugate`, `re`, `im`) can be
written out as CSV with the current implementation. But how will we do multiple
correlators? I see these options:

1.  We can just write out one CSV file per correlator matrix element. The
    information about $\vec d^2$, $\Gamma$, $\alpha$ and the $\{ \vec q_i \}$
    would then be encoded in the file name.

    This would be a classic thing, but these days I really dislike parsing
    filenames for information.

2.  Have one huge CSV file, where there are columns for everything. Then use
    *group by* and *summarize* operations to boil it down. The result would
    then also be a correlator in the *long data format*. This might be too
    inflexible to work with.

3.  Use a nested structure with a YAML representation. This way one can inject
    arbitrary amounts of meta data at every point in the structure. Also it
    could be parsed quite easily in R. But can it be generated from Wolfram
    Language? Apparently there are some third-party implementations, but the
    Wolfram Language can export JSON. I am perfectly fine with that, converting
    it to YAML is trivial with some Python or R script.

    The structure would look like this:
    \begin{multline*}
    \vec d \to \Gamma \to \alpha \to
    \Big\{
    \big(
    (O_\text{si}, O_\text{so},
    \vec q_\text{si,1}, \vec q_\text{si,2}, \ldots,
    \vec q_\text{so,1}, \vec q_\text{so,2}, \ldots),
    \\
    (\mathtt{datasetname}, \mathtt{conjugate}, \mathtt{re}, \mathtt{im})
    \big)
    \Big\} \,.
    \end{multline*}

As usual with these listings, I prefer the last option.

### Generating the structure

First we need to iterate through all the momenta $\vec d$ that we want to use.
For this we take the $\vec d_\text{ref}$ for each $\vec d^2$ and apply all
rotations from the quotient group $O_h / \mathrm{LG}(\vec d_\text{ref})$. This
way we get each unique $\vec d$ with the same magnitude. But in order to make
this easier, we just apply the whole of the octahedral group and just use
[`DeleteDuplicates`] afterwards.

This works fine for $\vec d^2 \leq 3$. With $\vec d = (0, 0, 2)$ we have the
same little group as for $\vec d = (0, 0, 1)$, but the $\vec d$ is different.
This needs to be addressed. As we are interested in low momenta at this point,
this corner is cut for now. We just use `Keys[IrrepDGammaAssoc[]]` to get all
the momenta that we have little groups for.

We just iterate through all $\vec d_\text{tot}$, then we iterate through all
irreps and then iterate through all relative momenta. Using our `MakeGroupSum`
function we construct the operators that will make up the correlator matrix.
Later on we need to take the outer product in order to obtain the actual
correlator matrix.

From the given operators we need to extract the momenta and then perform the
outer product of source and sink momenta. Then we need to apply the isospin
templates in order to get the actual HDF5 dataset names. From there we massage
the terms until we eventually get a JSON representation. This interface is
described in the next chapter.

---

-   **`UniqueTotalMomenta`**($\vec d^2$)

    Gives a list of unique total momenta corresponding to the given total
    momentum magnitude. For $\vec d^2 = 1, 2, 3, 4$, these are the following:

        {0, 0, 1}, {0, 0, -1}, {1, 0, 0}, {-1, 0, 0}, {0, 1, 0},
        {0, -1, 0}

        {1, 1, 0}, {1, -1, 0}, {-1, 1, 0}, {-1, -1, 0}, {0, 1, 1},
        {0, 1, -1}, {0, -1, 1}, {0, -1, -1}, {1, 0, 1}, {1, 0, -1},
        {-1, 0, -1}, {-1, 0, 1}

        {1, 1, 1}, {1, -1, -1}, {-1, 1, -1}, {-1, -1, 1}, {1, -1, 1},
        {1, 1, -1}, {-1, 1, 1}, {-1, -1, -1}

        {0, 0, 2}, {0, 0, -2}, {2, 0, 0}, {-2, 0, 0}, {0, 2, 0},
        {0, -2, 0}

-   **`RelativeToTotalMomenta`**($\vec d_\text{tot}$, $\{ \vec q_i \}$)

    Converts the given total momentum $\vec d_\text{tot}$ and $n - 1$ relative
    momenta $\{ \vec q_i \}$ to a list of $n$ particle momenta $\{ \vec p_i \}$
    using the relations
    $$ \vec p_1 = \vec d_\text{tot} - \sum_{j = 1}^{n - 1} \vec q_j
    \qquad \text{and} \qquad
    \vec p_i = \vec q_{i - 1} \,. $$

-   **`AllRelativeMomenta`**($\vec d_\text{tot}$, $\{ \vec q_i \}$, cutoff)

    If any of the momenta $\{ \vec p_i \}$ have a magnitude greater than
    *cutoff*, the result is an empty list.

-   **`MultiGroupSum`**(irrep, $\{\{ \vec p_i \}\}$, hold : [`Identity`])

    Calls `MakeGroupSum` for all the momentum sets given. Currently irrep row
    and column are fixed to $\alpha = \beta = 1$. Duplicates are automatically
    removed. Also this is limited to scalar particles with all $J_i = 0$ and
    $M_i = 0$.

    It can make sense to insert a [`Hold`] right before the `MakeGroupSum` such
    that you get a list of things that would be computed. As `MakeGroupSum`
    takes like 30 seconds per call, this function `MultiGroupSum` can take
    several hours to complete. If you pass [`Hold`] for *hold*, the actual spin
    projection will not be performed. Instead you will need to call `Map[expr,
    ReleaseHold, Infinity]` on the generated expression `expr` to evaluate all
    the `MakeGroupSum` calls.

-   **`GroupSumWholeIrrep`**($\vec d_\text{tot}$, irrep, $\{ \vec q_i \}$,
    cutoff, hold : [`Identity`])

    Performs the `MultiGroupSum` for all individual momenta that give the given
    total momentum.

-   **`GroupSumWholeTotalMomentum`**($\vec d_\text{tot}$, $\{ \vec q_i \}$,
    cutoff, hold : [`Identity`])

    Performs the `MultiGroupSum` for all irreps that are available with the
    given total momentum. Return value is an association from irrep to a list
    of spin projected operators.

# Design and implementation (R)

## Projecting the computed correlators

Taking the tables from the preceding step we must read in the prescribed HDF5
data sets and combine them given the factors. The interface that we get is the
following JSON structure:

```js
{
  "001": {
    "A1": {
      "000": {
        "000": {
          "C4cB_uuuu_g5.p00-1.d000_g5.p00-1.d000_g5.p000.d000_g5.p000.d000": {
            "conj": true, 
            "im": 0, 
            "re": 16
          }, 
```

It is an association of the following
$$ d_\text{tot} \to \Gamma \to O_i \to O_j \to C \to (\text{Re}, \text{Im}, \dagger) \,. $$

The levels are:

1. Total momentum (as string)
2. Irreducible representation
3. Label for the correlator matrix row, currently just a comma separated list
   of the relative momenta that have been used. In the future where we have
   other operators than just the pion, these will become more sophisticated
   labels. These are for human consumption anyway, so for the point of this
   code they are just strings.
4. Label for the correlator matrix column
5. HDF5 dataset name
6. Weight factor (`re` and `im`) as well as whether the numeric correlator
   needs to be complex conjugated first (`conj`)

### Iteration order

The contraction code generates one HDF5 file per diagram type and per gauge
configuration with names like `C4cC_cnfg4824.h5`. We have to decide how we
iterate through all this data. In the end we want one file per correlator
matrix element
containing all the observations from all configurations. These iteration orders
come to mind:

1.  Treat the configurations independent of each other. Open all the diagram
    files for a particular configuration (like 4825) and then build all the
    correlator matrix elements from this.

    *Advantages*: Only one HDF5 file has to be opened at any one time. Also the
    file can be consumed completely in one go, amortizing the expensive
    indexing operation after opening the file. Parallelization over
    observations is trivial this way.

    *Disadvantages*: The intermediate result are lots of correlator matrix
    elements but with only one observation each. These have to be combined
    later on. Depending on the data format it is just a concatenation. One has
    to be careful not to create millions of files for intermediate output.

2.  Focus on one correlator matrix element and build it for all observations at
    the same time.

    *Advantages*: There is no intermediate result, the whole statistics for a
    given correlator matrix element is created directly.

    *Disadvantages*: A lot of HDF5 files has to be opened, and only few
    datasets are going to be extracted. [We know
    that](https://github.com/HISKP-LQCD/hadron/issues/25) `h5ls` is about as
    expensive as `h5read`. Therefore this would make it rather costly to
    iterate this way. Also we would need a lot of RAM to store everything.

As discussed with Markus, the first option seems to be the more sensible one.
This decouples (numeric) projection and aggregation of the whole statistics. We
just need to think about the intermediate data format.

As we want to further process the data with R, there is no danger to use the
`Rdata` format for serialization. For each configuration the whole correlator
matrix will be
in exactly the same outer nested list structure, just with a numeric vector as
payload instead of the mapping between HDF5 dataset names and coefficients.

Combining these is simple: Just load them one at a time and do `rbind` on the
vectors. This will then be the data portion of an unsymmetrized `hadron::cf`
correlation function.

With this setup there will be no text IO for numeric data at all!

---

Tasks:

- Think about iteration order
- Iterate through files
- Open needed HDF5 files
- Combine various numeric correlators

## Comparison to Markus' data

Markus has already projected all the data for the $\rho$ channel, we can just
use this to compare. His data is not only in a different format but also uses
different parameterizations for individual momenta $p_i$ via his $P'$ and $q'$
values.

We want to do a fully automated complete comparison such that there are no
manual steps which are additional sources of errors.

### Momentum parameterization

Markus' parameterization works as
$$ p_1 = \frac{P'}{2} + q' \,, \quad p_2 = \frac{P'}{2} - q' \,. $$

The parameterization used in this code here is $p_1 = P - q$ and $p_2 = q$ for
the case of two particles. This means that in order to convert his
parameterization into mine one has to use the relations
$$ P = P' \,, \quad q = \frac{P'}{2} - q' \,. $$

The other way around works as
$$ P' = P \,, \quad q' = \frac{P - 2 q}{2} \,. $$

His correlator matrix format contains files called `operator-indices` that give
the three-vector $P$ that is actually used for that particular moving frame, as
well as the irrep row:

```
id	p_x	p_y	p_z	alpha
0	-1	0	0	1
1	0	-1	0	1
2	0	0	-1	1
3	0	0	1	1
4	0	1	0	1
5	1	0	0	1
```

One has to look up our $P$ and then find out the *operator id* from this table.

The correlator matrix elements are then labeled by their $q'_\text{source}$ and
$q'_\text{sink}$ values. These also correspond do indices that have to be read
off from the lines with the correct Dirac structure, `g: \gamma_{5},
\gamma_{5}`, in the file `gevp-indices`:

```
id	element
0	p: 1, g: \gamma_{50i}
1	p: 1, g: \gamma_{i}
2	p: 1, q: (0.0, 0.0, 0.5), g: \gamma_{5}, \gamma_{5}
3	p: 1, q: (0.0, 0.0, 1.5), g: \gamma_{5}, \gamma_{5}
4	p: 1, q: (1.0, 0.0, 0.5), g: \gamma_{5}, \gamma_{5}
5	p: 1, q: (1.0, 1.0, 0.5), g: \gamma_{5}, \gamma_{5}
```

From this we can deduce which correlator matrix rows and columns correspond to
my $q_\text{source}$ and $q_\text{sink}$.

### Momentum cutoff

For 000 T1u we are missing something in the HDF5 files:

```
Object C4cB_uuuu_p0-10.d000.g5_p-1-1-1.d000.g5_p111.d000.g5_p010.d000.g5
does not exist in this HDF5 file.
```

Also for 000 A2u we have this:

```
Object C4cB_uuuu_p-1-1-1.d000.g5_p-1-1-1.d000.g5_p111.d000.g5_p111.d000.g5
does not exist in this HDF5 file.
```

This is because Markus uses $d^2 \le 2$ for the rest frame. This is a some
additional rule that one has to keep in mind.

### Multiple irrep rows

With some d²=1 E we have this:

```
Error: nrow(filtered) == 1 is not TRUE
```

- 110
- 010
- 0-10
- 001
- 00-1

The issue here is that Markus has persisted multiple irrep rows. I need to
compare to all the ones that I have or just fix it to some particular one.

### Operators without coupling

Sometimes there are no operator indices:

```
cannot open file 'B35.32/rho_p3_A2_operator-indices.tsv': No such file or
directory
```

This likely is just because there is no coupling into that channel. My JSON
files also show that there are matrix elements but none of them contain any
actual correlators.

###

# Tests

The following are tests that we can perform to increase the confidence in $H_0$
which says that the code works perfectly fine.

## Physical non-coupling

There are a bunch of channels which should not couple with specific momenta and
irreps, we can check for these.

## Cross irrep

Correlation functions with operators from different irreps should just vanish.
We can just pick a few examples and see how that works out.

## Old projection code

Using Markus' projection code we can take correlators on a single gauge
configuration and project them to some irrep and momenta. The numeric results
should be exactly the same.

<!-- Links to Wolfram Language reference -->

[`Apply`]: http://reference.wolfram.com/language/ref/Apply.html
[`Association`]: http://reference.wolfram.com/language/ref/Association.html
[`ClebschGordan`]: http://reference.wolfram.com/language/ref/ClebschGordan.html
[`Coefficient`]: http://reference.wolfram.com/language/ref/Coefficient.html
[`Composition`]: http://reference.wolfram.com/language/ref/Composition.html
[`Dataset`]: http://reference.wolfram.com/language/ref/Dataset.html
[`DeleteDuplicates`]: http://reference.wolfram.com/language/ref/DeleteDuplicates.html
[`Dot`]: http://reference.wolfram.com/language/ref/Dot.html
[`EulerMatrix`]: http://reference.wolfram.com/language/ref/EulerMatrix.html
[`Function`]: http://reference.wolfram.com/language/ref/Function.html
[`Hold`]: http://reference.wolfram.com/language/ref/Hold.html
[`Identity`]: http://reference.wolfram.com/language/ref/Identity.html
[`List`]: http://reference.wolfram.com/language/ref/List.html
[`Map`]: http://reference.wolfram.com/language/ref/Map.html
[`NonCommutativeMultiply`]: http://reference.wolfram.com/language/ref/NonCommutativeMultiply.html
[`Part`]: http://reference.wolfram.com/language/ref/Part.html
[`Prefix`]: http://reference.wolfram.com/language/ref/Prefix.html
[`ReplaceAll`]: http://reference.wolfram.com/language/ref/ReplaceAll.html
[`Rule`]: http://reference.wolfram.com/language/ref/Rule.html
[`RuleDelayed`]: http://reference.wolfram.com/language/ref/RuleDelayed.html
[`StringJoin`]: http://reference.wolfram.com/language/ref/StringJoin.html
[`StringTemplate`]: http://reference.wolfram.com/language/ref/StringTemplate.html
[`TemplateApply`]: http://reference.wolfram.com/language/ref/TemplateApply.html
[`ToExpression`]: http://reference.wolfram.com/language/ref/ToExpression.html
[`WignerD`]: http://reference.wolfram.com/language/ref/WignerD.html
