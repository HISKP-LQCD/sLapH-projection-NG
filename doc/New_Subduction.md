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
  programming languages and libraries for the implementation. Also the design
  and implementation is covered here, including documentation for all functions
  in the Mathematica package.
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
  license riddled as well, at least *we* can use it.
  
    The octahedral group is not directly available in Mathematica but I managed
    to get it via the representation as a permutation group. I can also
    construct the little groups (also known as stabilizers) but not the irreps.

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
operators used and their correspondence in other languages. Mathematica
functions are in general hyperlinked directly to the official documentation.
You can also just type `?x` to get the help for `x`.

| Long | Infix | Python | R |
| --- | --- | --- | --- | --- |
| [`f`]`[x]` | `f @ x` or `x // f` | `f(x)` | `f(x)` |
| [`Apply`]`[f, a]` | `f @@ a` | `f(*a)` | `do.call(f, a)` |
| [`Map`]`[f, x]` | `f /@ x` | `map(f, x)` | `lapply(x, f)` |
| [`Composition`]`[f, g]` | `f @* g` | — | — |
| [`Function`]`[x, x^2]` | `#^2 &` | `lambda` | `function` |
| [`Part`]`[xs, i]` | `xs[[i]]` | `xs[i-1]` | `xs[[i]]` |
| [`Rule`] | `->` | — | — |
| [`RuleDelayed`] | `:>` | — | — |
| [`Association`]`[a -> b]` | `<| a -> b |>` | `{a: b}` | `list(a = b)` |
| [`Dot`]`[a, b]` | `a . b` | `a @ b` | `a %*% b` |
| [`NonCommutativeMultiply`] | `**` | — | — |
| [`ReplaceAll`] | `/.` | — | — |
| [`StringJoin`] | `<>` | `+` | 'paste0' |

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

Since the Wolfram Language is a functional one, the most straightforward way to
define a constant is as a function with zero arguments. In Haskell there is no
distinction between constants and argumentless functions at all.

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

    It seems that Mathematica has the notion of [`Dataset`], which means that
    one can just `Import` a CSV or TSV file and has a data frame as with R or
    Pandas. Also the `Dataset` supports having a complex matrix as column type,
    which would allow for a hybrid approach if we wanted to.

3.  There is another way, using nested [`Association`] instances. They are just
    like the `dict` in Python and the `std::map` in C++. This way we could make
    associations that map from some named group element to an association that
    maps from indices to a c-number.

    As associations can be called like functions, we can really nicely work
    with them. So let's make up a nonsense association like this:

    ```mathematica
    assoc = Association[foo -> bar, 1 -> "one", x^2 -> Association]
    ```

    When you now call `assoc[x^2]`, you get the function `Association`. Or you
    call it with `assoc[1]` and get `"one"`. Missing keys are some sort of
    error condition via say `Missing["KeyAbsent", bar]`.

Markus was immediately sold on this approach 3, and I quite like it as well.
This is now implemented.

---

-   **`ReadDataframe`**(filename)

    Reads a CSV formatted table with a header and generates a `Dataset` from
    it.

-   **`ConvertMatrixElement`**(dataset)

    The matrix elements of the irreducible representations,
    $D^\Gamma_{\alpha\beta}(g)$ are stored as strings in the CSV file. Using
    `ToExpression` these are parsed into Mathematica symbols. This function
    takes the dataset and replaces the “matrix element” column with the parsed
    one.

-   **`ReadIrreps`**(filename)

    The composition of `ConvertMatrixElement` and `ReadDataframe`.

-   **`MatrixElementsToAssociation`**(dataset)

    Converts the “row”, “col” and “matrix element” into mappings of the form
    $(\alpha, \beta) \to D^\Gamma_{\alpha\beta}(g)$. The result is a `Dataset`
    with just one row per irrep and little group element and an `Association`
    which contains all the matrix elements.

-   **`IrrepsToAssociation`**(matrixElementsAssociations)

    Continuing with the result from `MatrixElementsToAssociation` this function
    then builds a more deeply nested `Association` which has the form
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
    momentum $\vec p_\text{cm}$. As this function is implmented via
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
Exchange](https://physics.stackexchange.com/questions/459408/clebsch-gordan-coefficients-for-more-than-2-particles), and in one answer I was told that this is correct though not unique.

I was pointed to the [Racah $W$
coefficient](https://en.wikipedia.org/wiki/Racah_W-coefficient). These do not
seem to be implemented in Mathematica, though.

For coupling pions these factors will all be 1 anyway, so we could consider not
bothering with them for just now.

---

-   **`HigherClebschGordan`**($\{ j_i \}$, $\{ m_i \}$, $(J, M)$)

    Computes the Clebsch-Gordan coefficient for coupling all the spins $(j_i,
    m_i)$ together to $(J, M)$. Note that the API is different from
    `ClebschGordan` as here the $j_i$ and $m_i$ are grouped among each other
    instead of grouped by $i$.

    In contrast to `ClebschGordan` this function does not emit any warnings in
    case the coupling is unphysical.

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

-   **`ExtractMomenta`**(expr)

    Extracts all the momenta from the `SingleOperator`s used. The resulting
    expression contains summand of this form:

    ```mathematica```
    (-(1/2) + (I Sqrt[3])/2) {1, 0, 1} ** {0, 1, 0}
    ```

    By replacing the [`NonCommutativeMultiply`] with something else, one can
    extract just the momenta from this. One has to be a bit careful as
    Mathematica will try to distribute the common factor onto the elements in
    the lists if one replaces the [`NonCommutativeMultiply`] with something
    simple like [`List`].

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

Similarly we define a $\pi^-$ operator:

```mathematica
\[Pi]Minus[s1_, s2_, c1_, x1_] :=
    FieldB[dn, c1, s1, x1] ** (Gamma^5)[SI[{s1, s2}]] ** 
    Field[up, c1, s2, x1];
```

From these we build an $I = 2$ operator.

```mathematica
\[Pi]\[Pi]I2[s1_, s2_, s3_, s4_, c1_, c2_, x1_, x2_] :=
    \[Pi]Plus[s1, s2, c1, x1] ** \[Pi]Plus[s3, s4, c2, x2]; 
```

For reasons not entirely clear to me we need to define the “bar” version of
that operator ourselves. Perhaps it is because the “quark contraction tool”
does not know about the flavors? Anyway, the second operator is this:

```mathematica
\[Pi]\[Pi]I2Bar[s1_, s2_, s3_, s4_, c1_, c2_, x1_, x2_] :=
    \[Pi]Minus[s3, s4, c2, x2] ** \[Pi]Minus[s1, s2, c1, x1]; 
```

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
that looks like this:
\begin{align*}
&\text{tr}\left(\gamma^5\;S^{\text{up}}(p_4,p_1)\;\gamma^5\;S^{\text{dn}}( p_1,p_4)\right) \text{tr}\left(\gamma^5\;S^{\text{up}}(p_3,p_2)\;\gamma^5\;S^{\text{dn }}(p_2,p_3)\right)
\\&\quad
+\text{tr}\left(\gamma^5\;S^{\text{up}}(p_3,p_1)\;\gamma^5\;S^{\text{dn}}(p_1,p_3)\right) \text{tr}\left(\gamma^5\;S^{\text{up}}(p_4,p_2)\;\gamma^5\;S^{\text{dn }}(p_2,p_4)\right)
\\&\quad
-\text{tr}\left(\gamma^5\;S^{\text{dn}}(p_2,p_3)\;\gamma^5\;S^{\text{up}}(p_4,p_2)\;\gamma^5\;S^{\text{dn}}(p_1, p_4)\;\gamma^5\;S^{\text{up}}(p_3,p_1)\right)
\\&\quad
-\text{tr}\left(\gamma^5\;S^{\text{up}}(p_4,p_1)\;\gamma^5\;S^{\text{dn}}(p_2,p_4)\;\gamma^5\;S^{\text{up}}(p_3,p_2)\;\gamma^5\;S^{\text{dn}}(p_1,p_3)\right)
\end{align*}

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

The trace is cyclic. Therefore similar expressions can be equivalent. The trace
expression used, `trace` from `qct`, does not let Mathematica know that it is
cyclic. Therefore we manually cycle them until the propagator starting from the
lowest vertex (to be defined) is the second element. The first element shall be
the $\Gamma$-structure. The first vertex is the source vertex with the lowest
number. If there is no source vertex, then the lowest sink vertex is to be
used.

The key observation is that [expressions can be used like
lists](https://reference.wolfram.com/language/tutorial/PartsOfExpressions.html).
With this we can access the different parts of expressions. The arguments of a
function are subexpressions.

---

-   **`RotateGammaToFront`**(expression)

    Takes an expression like 

    ```mathematica
    Gamma^5.DE[{"dn", "dn"}, {so[2], si[1]}].
        Gamma^5.DE[{"up", "up"}, {si[2], so[2]}]
    ```

    and cyclicly permutes the expression if it does not start with a `Gamma`.
    The resulting expression always start with `Gamma`.

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

### Contractions into dataset name templates

At the stage of the Wick contractions we still know which labels correspond to
source and sink vertices. We can track this information by adding the undefined
functions `so` and `si` which label the source and sink vertices instead of
`x1` and so on. The contractions are then defined with operators evaluated like
so:

```mathematica
wc = WickContract[\[Pi]\[Pi]I2Bar[s1, s2, s3, s4, c1, c2, so[1], so[2]] **
    \[Pi]\[Pi]I2[s5, s6, s7, s8, c5, c6, si[1], si[2]]];
```

In the contracted form we have propagator expressions like these:

```mathematica
DE[{up, up}, {si[1], so[2]}
```

These have to be mapped to HDF5 dataset names. Actually they have to be mapped
to templates ([`StringTemplate`]) that can then be used with [`TemplateApply`]
to insert the momenta. The following replacement can be used to replace an
expression with a `C4cD` diagram:


```mathematica
qc /. trace[Gamma^g1_ .DE[{f1_, f1_}, {si[si1_], so[so2_]}].
     Gamma^g2_ .DE[{f2_, f2_}, {so[so2_], si[si1_]}]] trace[
    Gamma^g3_ .DE[{f1_, f1_}, {si[si3_], so[so4_]}].
     Gamma^g4_ .DE[{f2_, f2_}, {so[so4_], si[si3_]}]] :> 
  TemplateApply[
   "C4cD_uuuu_g`g1`.p`x1`.d000_g`g2`.p`x2`.d000_g`g3`.p`x3`.d000_g`g4`\
        .p`x4`.d000",
    <|"g1" -> g1, "x1" -> "`pso" <> ToString @ so2 <> "`",
      "g2" -> g2, "x2" -> "`psi" <> ToString @ si1 <> "`",
      "g3" -> g3, "x3" -> "`pso" <> ToString @ so4 <> "`",
      "g4" -> g4, "x4" -> "`psi" <> ToString @ si3 <> "`"|>]
```

This way the we replace the two summands

```mathematica
trace[Gamma^5.DE[{up, up}, {si[1], so[2]}].
    Gamma^5.DE[{dn, dn}, {so[2], si[1]}]]
  trace[ Gamma^5.DE[{up, up}, {si[2], so[1]}]..
    Gamma^5.DE[{dn, dn}, {so[1], si[2]}]] + 
trace[Gamma^5.DE[{up, up}, {si[1], so[1]}].
    Gamma^5.DE[{dn, dn}, {so[1], si[1]}]]
  trace[Gamma^5.DE[{up, up}, {si[2], so[2]}].
    Gamma^5.DE[{dn, dn}, {so[2], si[2]}]]
```

with this:

```mathematica
"C4cD_uuuu_g5.p`pso1`.d000_g5.p`psi1`.d000_g5.p`pso2`.d000_g5.p`psi2`.d000" +
  "C4cD_uuuu_g5.p`pso2`.d000_g5.p`psi1`.d000_g5.p`pso1`.d000_g5.p`psi2`.d000"
```

It is a sum of two strings, which is exactly what we want. In a later stage we
can insert the momenta corresponding to the spin projected operators.

---

-   **`MakeTemplate`**($n$)

    Creates a template for an $n$-particle correlation function. For the case
    of $n = 2$, the result looks like this:

    ```mathematica
    "g`g1`.p`x1`.d000_g`g2`.p`x2`.d000"
    ```

-   **`DatasetNameRules`**()

    List of delayed rules ([`RuleDelayed`]) that can be used with with
    [`ReplaceAll`] to replace propagator expression with HDF5 dataset name
    templates ([`StringTemplate`]). The replacements are expressions like
    these:
    
    ```mathematica
    "C2cD_uu_g5.p`pso1`.d000_g5.p`psi1`.d000"
    ```

## Wick contraction and spin

We now have an operator that is projected with spin. And we have the Wick
contraction which already has the isospin incorporated into it. One big
difference is that the spin operator is at sink (or source) only, whereas the
Wick contraction is a full correlation function already.

As a next step we just build a full correlation function from the spin
operators, namely writing $O O^\dagger$. This has all the factors that we need,
we just have to replace the formed products of single-particle creation and
annihilation operators with the output from the Wick contraction.

We will use pattern replacement for that. As a toy example, we take this
expression:

```mathematica
toy = f[x] ** f[y];
```

The `f` stands for the function `SingleOperator` with all its arguments. We
want to replace it with the function `g`, which stands for the multi-particle
operator. We can do this directly with a pattern replacement:

```mathematica
toy /. f[a_] ** f[b_] -> g[a, b]
```

The result is `g[x, y]`, as desired. The same works with more complicated
expressions with more factors.

The replacement that we want to apply is the following:

```mathematica
repl :=  qc /. x1 -> #1 /. x2 -> #2 /. x3 -> #3 /. x4 -> #4 &
```

This is a function that takes four momenta. It then takes the quark contraction
`qc` and replace the four position/momentum labels with the arguments. The
function `ReplaceSingleOperatorScalarWick` accepts the correlator of spin
operators and this replacement function and gives us the variant with concrete
momenta. We obtain
\begin{align*}
&
\text{tr}\left(\gamma^5\;S^{\text{up}}(\{0,1,1\},\{1,0,0\})\;\gamma^5\;S^{\text{dn}}(\{1,0,0\},\{0,1,1\})\right)
\\&\qquad\times
\text{tr}\left(\gamma^5\;S^{\text{up}}(\{1,0,0\},\{0,1,1\})\;\gamma^5\;S^{\text{dn}}(\{0,1,1\},\{1,0,0\})\right)
\\&\quad
+\text{tr}\left(\gamma^5\;S^{\text{up}}(\{0,1,1\},\{0,1,1\})\;\gamma^5\;S^{\text{dn}}(\{0,1,1\},\{0,1,1\})\right)
\\&\qquad\times
\text{tr}\left(\gamma^5\;S^{\text{up}}(\{1,0,0\},\{1,0,0\})\;\gamma^5\;S^{\text{dn}}(\{1,0,0\},\{1,0,0\})\right)
\\&\quad
-\text{tr}\big(\gamma^5\;S^{\text{dn}}(\{1,0,0\},\{0,1,1\})\;\gamma^5\;S^{\text{up}}(\{1,0,0\},\{1,0,0\})
\\&\qquad\times
\gamma^5\;S^{\text{dn}}(\{0,1,1\},\{1,0,0\})\;\gamma^5\;S^{\text{up}}(\{0,1,1\},\{0,1,1\})\big)
\\&\quad
-\text{tr}\big(\gamma^5\;S^{\text{up}}(\{1,0,0\},\{0,1,1\})\;\gamma^5\;S^{\text{dn}}(\{1,0,0\},\{1,0,0\})
\\&\qquad\times
\gamma^5\;S^{\text{up}}(\{0,1,1\},\{1,0,0\})\;\gamma^5\;S^{\text{dn}}(\{0,1,1\},\{0,1,1\})\big)
\,.
\end{align*}

This is now very close to the stuff that we actually compute, we should be able
to get HDF5 dataset names out of that expression.

---

-   **`ReplaceSingleOperatorScalar`**(expression, replacement)

    Takes the *expression* stemming from `MakeGroupSum` for just scalar
    particles and replaces all the occurences of
    `ConjugateTranspose[SingleOperator[…]]` and replaces it with a single
    multi-particle operator. The *replacement* has to be a function taking the
    two momentum arguments.

    For the $\rho$-channel we use a two-pion operator and therefore our
    *replacement* is the following function:

    ```mathematica
    \[Pi]\[Pi]I1[s1, s2, s3, s4, c1, c2, #1, #2] &
    ```

-   **`ReplaceSingleOperatorScalarWick`**(expression, replacement)

    Replaces *expression* from two to two particles like this one

    ```mathematica
    ConjugateTranspose[SingleOperator[1, 0, 0, p1_]] ** 
    ConjugateTranspose[SingleOperator[2, 0, 0, p2_]] ** 
    ConjugateTranspose[
      ConjugateTranspose[SingleOperator[1, 0, 0, p3_]] **
      ConjugateTranspose[SingleOperator[2, 0, 0, p4_]]]
    ```

    with the *replacement* being evaluated with the momenta as:

    ```
    replacement[p1, p2, p3, p4]
    ```

-   **`MomentumToString`**($\vec p$)

    Converts an integer momentum vector like `{1, 0, -1}` into the string
    `"10-1"`.

-   **`MomentaToRules`**(momenta, location)

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

-   **`MomentumPluginRecursive`**(rules, templateExpression)

    Recursively traverses through the *template expression* and replaces the
    momenta with the rules. The *rules* must be an [`Association`] of this form:

    ```mathematica
    <|"psi1" -> "011", "psi2" -> "100",
      "pso1" -> "011", "pso2" -> "100"|>
    ```

# Tests

- Physical non-coupling
- Cross-Irrep

<!-- Links to Wolfram Language reference -->

[`Apply`]: http://reference.wolfram.com/language/ref/Apply.html
[`Association`]: http://reference.wolfram.com/language/ref/Association.html
[`ClebschGordan`]: http://reference.wolfram.com/language/ref/ClebschGordan.html
[`Composition`]: http://reference.wolfram.com/language/ref/Composition.html
[`Dataset`]: http://reference.wolfram.com/language/ref/Dataset.html
[`Dot`]: http://reference.wolfram.com/language/ref/Dot.html
[`EulerMatrix`]: http://reference.wolfram.com/language/ref/EulerMatrix.html
[`Function`]: http://reference.wolfram.com/language/ref/Function.html
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
[`WignerD`]: http://reference.wolfram.com/language/ref/WignerD.html
