% Nuweb formatted latex file 
% Most of this is bog standard latex with fortran or c rolled in
% Anything to do with @@ characters is probably specific to nuweb
%
%
% The word FIXME anywhere in this document indicates 
% an area where more attention is still needed.
%
%  Copyright (C) 2004 by Jonathan Wright
%     Grenoble, France
%     email: wright@@esrf.fr
%
%   This file is part of Linarp.
%
%   Linarp is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation; either version 2 of the License, or
%   (at your option) any later version.
%
%   Linarp is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with Linarp; if not, write to the Free Software
%   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%
\documentclass[10pt,a4paper,twoside,notitlepage]{report}

\usepackage{graphics} % For the pictures
\usepackage{anysize}  % Try to circumvent Latex default margins
\usepackage{fancyhdr}
\usepackage[dvipdfm,bookmarks=true,backref,bookmarksnumbered=true,bookmarkstype=toc]{hyperref}
\newcommand{\var}[1]{\textbf{\textsf{#1}}} % highlight variables in text
\newcommand{\code}[1]{\textbf{\textsf{#1}}} % highlight code in text
\newcommand{\param}[1]{\textbf{\textsf{#1}}} % ... parameters ...
\newcommand{\mod}[1]{\textbf{\textsf{#1}}} % ... module names ...
\newcommand{\mb}  [1] {\mathbf{#1}}


\begin{document}

%\marginsize{1.5cm}{1.5cm}{1.5cm}{1.5cm} % Needs anysize
%\pagestyle{headings}            % These are ugly - fix them somehow?

\pagestyle{fancy}
\renewcommand{\chaptermark}[1]{
      \markboth{\chaptername
      \ \thechapter.\ #1} {} }

\renewcommand{\sectionmark}[1]{
      \markright {   
      \ \thesection.\ #1} {} }

\fancyhead[LE,RO]{\rightmark}
\fancyhead[LO,RE]{\leftmark}
\fancyfoot[C]{\today}
\fancyfoot[LE,RO]{\thepage}
\fancyfoot[LO,RE]{J. P. Wright}
\renewcommand{\footrulewidth}{0.4pt}

\pagenumbering{arabic}          % Page numbers



\title{\textbf{\textsf{Linarp}} \\ Linarp Is Not A Rietveld Program}
\author{JPW + anyone who cares to contribute}
\date{Started May 4, 2004, already it is \today}

\maketitle

\abstract{A two step approach to crystal structure refinement from powder
diffraction data is of interest to the author. 
This program aims to carry out refinements of crystallographic models
where the integration of the powder data and the structure refinement
are separated.
The name Linarp derives from the recursive acronym that 
\textbf{L}inarp \textbf{I}s \textbf{N}ot \textbf{A} \textbf{R}ietveld 
\textbf{P}rogram
}

\tableofcontents

FIXME

\begin{itemize}
\item add a preface to explain about the document, how to read
        it and modify the program
\item What does it mean when polynomials are orthogonal?
\item Set up the "unittest" python module for code testing
\item Set up the "doctest" for ditto
\item Add the "testdenseleastsquares.py" to automated testing framework
\item Look into the inspect module for wrapping a gui to arbirary objects
\item Make the testdenseleastsquares a class which is the basis for refinements?

\item Split the peak fitting python thing into a 1D data weight matrix and 
      dense gradients optimiser
\item Make a dense gradients function thing from the peak fitter.
\item Make an analogous 2D weight matrix optimiser for the bulk solvent
      refinement and apply it to a protein and LaB6
\item Make a sparse gradients 1D data object for pattern decompositions
\item Make the gui be constructed from a configuration file to have it
      easier to upgrade/change/modify etc.
\end{itemize}

\section*{Index of file names}
@f 
\section*{Index of macro names}
@m 

\section*{Preamble}

If you are reading this document, and you did not write it, then quite
a lot of it will not make much sense.
These are the rambings of a slightly deranged crystallographer, who has 
been thinking out loud.
Hopefully the text will be tidied up before any of this is released
to the general public, but if you have gotten hold of a copy in the meantime,
then please accept my apologies for any pain incurred when reading ahead.
The idea is just to make a refinement package which does everything
I can think of quickly and correctly, and will be easy to modify
in the future for the things I haven't thought of yet.

The text is split into chapters. 
The introduction gives an overview for the main ideas which motivate
this work.
There is then a chapter which attempts to build a computing
system from the various components - aiming to define interfaces
between the different parts of the code. 
We don't want to have to rewrite any part of the code when we change
our mind about some area of it.
Then there are specialised chapters for each area of the code, things
on peakshapes, crystal structures, symmetry, etc.
Eventually there should be a chapter on the gui part of the code, so
that the methods developed can interact with humans too.



\chapter{Introduction}

JPW very recently wrote a paper which drew heavily on previous work
by others which showed that the two step approach to refinements of
models from powder data could offer some advantages. 
Actually carrying out such refinements will need a refinement program.
The sort of thing that other people can use too.
Making such a program is an extremely ambitious goal, particularly for
one person who expects the program to do absolutely everything.
Nevertheless, by continuing to draw heavily on the work of others,
it might be possible to get something up and running in a reasonable
timescale.

Eventually it might be nice to have a book which describes exactly
how the program works, documenting the implementations and
mathematics needed.
In many ways the documentation would be of more value than the program;
in that generating the book is far more educational and more likely to 
lead to a program which is correct.

\section{Goals}

The first milestone - a refinement using a non-diagonal weight matrix
has already been achieved. 
Unfortunately the implementation of a hair raising python script
to carry out a refinement is not going to be acceptable for anyone,
including the author, in the future.
A more general program is needed. 
Some of the things the program should do would be the following:

\begin{itemize}
\item The 3 ``M''s
 \begin{itemize}
 \item Multiple crystallographic phases
 \item Multiple histograms of data
 \item Magnetism
 \end{itemize}
\item Correct propagation of error bars onto derived parameters
\item Generation of arbitrary constraints and restraints
 (eg NCS, adp's to respect bondlengths)
\item Sanity preserving optimisation procedures
\item The ability to correctly apply arbitrary symmetry 
(from generators) to arbitrary properties
\item A modular, extensible design which does not 
seriously compromise efficiency
\end{itemize}

\section{Realistic plans}

Most interesting at the moment is to use the non-diagonal matrix
approach for refinement of protein structures from powder data.
A first step would be to actually refine protein structures from
powder data via this approach.
In order to get the ambitious project to eventually lift off
of the ground, it would be handy to make some key design decisions 
(bear in mind that these will be scrapped and the whole thing rewritten
from scratch at a later date).

The next measurable goal is therefore to optimise a protein
structure against powder data, using some system of restraints.
This goal is set bearing in mind that the MMTK package for
python already optimises protein structures without any data
present, so it is just a question of supplying derivatives
with respect to the powder data to that package.
Structure factors, derivatives and symmetry are all 
available in the cctbx package.
With the things which are already available, this seems
like a goal which can be obtained on a reasonable timescale.
A few days for a one off, maybe a little longer for something
which can grow in the future.

Later milestones should include a more user friendly method 
for obtaining correlated integrated intensities. 
Currently the PRODD computer program can carry out this task,
and has already been interfaced to python via the writing out
of a suitable file.

\section{Mathematical background}

In thinking about the design of the program it is going
to be worthwhile to run through the basics of the mathematics
that we intend to implement.
The fundamental idea is to transform the loops inside a Rietveld
program so that the structure refinement can be completely
separated from the pattern fitting.

To do that you need to know how a Rietveld program 
(is supposed to) work in the 
first place.
This is how one program (PRODD) works, I guess others are 
probably similar, based on their behaviour.
You read in a structure and dataset. 
Then generate a calculated curve based on the structure, 
and derivatives of that curve at each datapoint.
The model is:
\[ y_{c,i} = b_i(\mb{x}) +
    \sum_\mb{H} I_\mb{H}(\mb{x}) \phi_{i,\mb{H}}(\mb{x}),\]
where $\mb{x}$ is the vector of parameters to optimise. 
$b$ is the backgrouns, $I$ is the peak intensity and $\phi$
picks up the peakshape, scale factor and any other corrections
needed to turn a peak into diffractometer counts. $\mb{H}$
is the hkl of the reflection, and the summation is over contributing
peaks to datapoint $i$ (or $\phi$ is zero for
those that don't contribute and your program is very slow).
The least squares problem is then:
\[ \chi^{2}_{R} = \sum_i \frac{ \left( y_{o,i} - y_{c,i}(\mb{x}) \right) ^2 }  
{\sigma_i^2} \]
Our model is linearised by Taylor expansion to first derivatives.
\[  y_{c,i}(\mb{x}) =  y_{c,i}(\mb{a}) + 
              \sum_k s_k \frac{\partial  y_{c,i}(\mb{a})}{\partial x_k} \]
An iterative algorithm is then used which solves the series
of linear problems which are given by adding the shifts 
($s_k=x_k-a_k$).
Defining $ d_i = y_{o,i} - y_{c,i}(\mb{x}) $ and 
$ g_{ik} = \partial  y_{c,i}(\mb{a})/\partial x_k $
and substituting the
Taylor expansion into the $\chi^2_R$ function gives:
\[ \chi^{2}_{R} = \sum_i \frac{ \left(
d_i -    \sum_k s_k g_{ik}
 \right) ^2 } {\sigma_i^2} \]
Differentiating with respect to shift $s_m$ is going to 
give us a least squares problem to solve:
\[ \frac{\partial \chi^2_R}{\partial s_m} =
   \sum_i -2 \frac{g_{im}} {\sigma_i^2} 
   \left(d_i - \sum_k s_k g_{ik} \right) =0 \]
At a minimum the derivatives with respect to any parameter are
zero and so we get:
\[ \sum_i \frac{g_{ik} d_i}{\sigma^2_i} = 
   \sum_{ik} \frac{g_{im} g_{ik}}{\sigma^2_i} s_k \]
(both sides were divided by 2).
There are as many equations as there are parameters to find and
we have the normal equations to solve.
Identifying $[\beta]_m = \sum_i d_ig_{im}/\sigma^2_i $ and 
$[\alpha]_{km} = \sum_{i} g_{im}g_{ik}/\sigma^2_i $
is going to allow us to write that in compact notation as:
\[ [\alpha]\mb{s} = [\beta] \]
Convince yourself. 
The interesting part is the generation of the least squares
matrix elements $ [ \alpha ]_{km} $. 

In the ``conventional'' Rietveld approach we compute the
actual values of these numbers by summing over datapoints.
Thus:
\[ [\alpha]_{km} = \sum_i \frac{1}{\sigma^2_i}
          \frac{\partial  y_{c,i}(\mb{a})}{\partial s_k}
          \frac{\partial  y_{c,i}(\mb{a})}{\partial s_m} \]
Going back to our model for the data, the derivative
terms are actually coming from a sum of contributions to
the datapoint (some peaks and a background).
Generalising beyond peak intensities and background to 
call it a sum of linear functions, $f$, the model is then:
\[ y_{c,i}(\mb{a}) = \sum_p f_p(\mb{a}). \]
This makes background be treated as just another peak 
to add to the summation.
The Taylor expansion gives:
\[ y_{c,i}(\mb{x}) = \sum_p f_p(\mb{a}) + 
 \sum_{pk} s_k \frac{\partial f_p(\mb{a})}{ \partial s_k} \]
For the peaks in the model, the functions $f_p$ are given by:
\[ f_p(\mb{a}) = I_\mb{H}(\mb{a}) \phi_\mb{H}(\mb{a}) \]
So they are the product of two things. The actual Bragg
peak intensity as you might manage to measure with a single
crystal, and a peak shape function $\phi$ which converts
an overall peak intensity into a number of diffractometer 
counts, taking into account all the polarisation, lorentz
and other effects which prevent humans doing crystallography
with pencils and paper.
Derivatives of $f_p$ with respect to any parameters $s_k$
are going to be given by:
\[ \frac{\partial f_p(\mb{a})}{\partial s_k} =
  I_\mb{H}(\mb{a})  \frac{\partial \phi_{\mb{H}}(\mb{a})}{\partial s_k} +
  \phi_\mb{H}(\mb{a}) \frac{\partial I_{\mb{H}}(\mb{a})}{\partial s_k} \]
At this point, in practice, the comments in at least one code
claim to be using the chain rule here, as they have
noticed that one of the two terms is always zero for most
commonly occuring derivatives. (Is there something where non-zero
derivatives are needed for both?).
For example, if you refine the positional parameters for a crystal
structure, then the peak shape function does not change.
The trick, is to try to partition the model into two sets
according to which of the terms in this equation is the non-zero
one. Only cell parameters (via distance restraints) 
are really trying to shoot us in the foot, and only protein
crystallographers currently seem to have noticed this (there are
validation packages which refine the cell parameters to fit restraints
better when you have high resolution data and the atomic parameters
were not restrained in the structure refinements).

So getting back on track, we had $ [\alpha]_{km} $ as a sum over data
points. 
In practice it is actually a sum over contributions to data points
and the purpose of this whole program is just to interchange
the order of summations:
\[ [\alpha]_{km} = \sum_i  \frac{1}{\sigma^2_i} 
  \left( \sum_{p} \frac{\partial f_p(\mb{a})}{ \partial s_k} \right)
  \left( \sum_{q} \frac{\partial f_q(\mb{a})}{ \partial s_k} \right)
 \]
Considering only structural parameters gives:
\[ [\alpha]_{km} = \sum_i  \frac{1}{\sigma^2_i} 
  \left( \sum_{p} \phi_{p,i}\frac{\partial I_p(\mb{a})}{ \partial s_k} \right)
  \left( \sum_{q} \phi_{q,i}\frac{\partial I_q(\mb{a})}{ \partial s_k} \right)
\]
\[ [\alpha]_{km} = \sum_{pq} 
  \frac{\partial I_p(\mb{a})}{ \partial s_k} 
  \frac{\partial I_q(\mb{a})}{ \partial s_k} 
  \sum_i \frac{\phi_{p,i}\phi_{q,i}}  {\sigma^2_i}  
\]
So we can compute the peak shape function once, and fill out an
npeaks $\times$ npeaks matrix of values in that summation over $i$.
That matrix has a special sparse structure, in that only overlapping
contributions have non zero matrix elements
($\phi$ is zero a lot of the time).
The idea is then to record this ``non-diagonal-weight-matrix'' so
as to avoid a lot of tedious calculations in trying out structural models.
Also we want to be routinely comparing the structure we have to the one
which has the maximum possible number of free parameters (all of the peak
intensities as variables). 
How well do we do compared to a pattern composition gives a better idea
about the quality of the structure compared to how well do we do in
an absolute sense.
The rhs of least squares problem comes from the residuals in the current
model. 
(That's a detail - you can also go directly for a correlated $\chi^2$ 
function and use $I_{obs}-I_{calc}$ where $I_obs$ comes from the
pattern decomposition).

\section{Partitioning models}

The main goal of the project is to \emph{strictly}
partition a refinement into the things which are to
do with intensity and the things which are not.
For any point in the raw data we have:
\[ \frac{\partial y_c}{\partial x} = \sum
  \phi \frac{\partial I}{\partial x} +
  I \frac{\partial \phi}{\partial x} \]








\section{Design considerations}

Python has become the author's programming language of choice. 
It allows both c and fortran to be used together in some
sort of harmony and gives a nice way of throwing things
together in a very high level language.
Several things which need to done are forseen at the moment,
some of these would ideally share the same code base, hence
the enumeration of things to do, with some thought as to 
how to do them.

\begin{itemize}
\item Fit the diffraction pattern(s) using up-to-date peakshapes
\item Fit models to the CII data
\item Do clever things - like maximising entropy, likelihood etc
\item Produce fourier and patterson maps
\item Display structures, restraints and maps
\item Derive correctly symmetry constraints
\item Determine esd's on parameters and derived parameters
\end{itemize}

That list is not exhaustive, but it should give some overall ideas
of what is going to be eventually needed.
The verb ``fit'' implies that the least-squares come optimisation 
module would ideally be shared between a pattern fitting program
and a structure optimisation program.
This would want to allow arbitrary things like ``entropy'' to be
taken into account somehow.
Display means getting to grips with some graphics - 2D for plotting
fits to data and 3D for displaying structures (matplotlib and
pymol?).

In order for all of this to work, it seems there will be a need
to pass numbers between different algorithms which can 
perform calculations with them.
Given that python is used there is a temptation to use
some modern data structures and then write complex wrappers
to be able to access these from c/fortran.
Alternatively we can use c/fortran data structures (arrays as raw 
strings) and let python scratch it's head to read them.
From the efficiency point of view, we would like python to never
do much more than supply pointers to extension routines, who will
know how to access and process any data.
The Numeric (or numarray) packages would seem like a good starting place
for passing chunks of numbers around, provided we have a way of 
ensuring that the type and memory layout are the same as that expected
by the c/fortran routines.
So, our basic data structure is going to a Numeric (or numarray) array
and the dimensions and type are going to be decided so that they can
be passed into fortran routines without generating temporaries.

What are the data structures going to need to represent and how should
they represent them?
Some things to represent are:
\begin{itemize}
\item Crystal structures
 \begin{itemize}
 \item Atomic positions and properties (type or scattering factor, magnetic moment,
  thermal pars etc)
 \item Symmetry operators
 \item Reflections (hkl indices, structure factors)
 \end{itemize}
\item Data
 \begin{itemize}
 \item Raw experimental data - counts and experiment geometry
 \item Statistic weights
 \item Model for the data  - peak shapes and positions
 \end{itemize}
\item{Least squares problems}
 \begin{itemize}
 \item Models to compute the data
 \item Pseudo-data to restrain the models
 \item Parameters to be optimised
 \item Relationships between parameters in the model and parameters in the 
       optimisation problem (constraints)
 \item Statistical weights on the data
 \item Error bars for the refined parameters
 \item Methods for propagating (correctly) the errors on the refined
       parameters into errors on any derived parameters
 \end{itemize}
\end{itemize}

Some of these things should obviously be represented by simplistic arrays.
For example hkl indices and raw data fit quite nicely into arrays.
Symmetry operators can also be represented as matrices ($3\times3$~for rotations
plus a translation vector. Some trick exists with $4\times4$~arrays for doing
that?
Properties like magnetic moments and thermal factors depend on the symmetry
of the atomic site. 
Peak shapes depend on hkl, the symmetry and the experimental geometry, but I guess
not the atomic positions.
For magnetism and incommensurate structures the hkl lists need more than 3 indices.
The least squares optimiser needs to talk to these various components 
in an intelligent way.

Despite that Linarp is not intended to be a Rietveld program, it ought
to be able to do Rietveld refinement in a transparent way. 
This is just a question of doing some summations in a different order.

\chapter{Data structures and algorithms}

Before getting something strung together and working, I want to decide how
the various bits of the program are going to communicate with each other.
If they can all just access a particular set of data structures then life
will be much easier.
Experience has shown that most of my difficulties with things like this have
been in passing around data between different bits of programs.
Scaling has certainly suffered from poor initial choices.

The pattern fitting will need to know the hkl lists for each of the 
crystal structures which are present. 
Logically it would want them to be sorted into order of d-spacings 
with different crysalline phases merged together to give a single list.
So far as the crystal structure is concerned, it only needs the d-spacings,
hkl's and crystal symmetry and unit cell. 
There is no need to know anything about atomic positions.
Of course, all of the instrumental parameters are needed here.
This program should be able to fit normal one dimensional datasets, but
it would be good to make it general enough to fit higher dimensional
datasets as well, without changing the underlying implementation.
Inspired by CCSL - this would be a peak centre function and a peak shape
function.
The pattern fitting may optionally be supplied with a list of intensities
to use in the fitting.

Crystal structure fitting needs to know the about the kind of radiation
used to get scattering factors.
It also needs to know about the weight matrix and extracted intensities.
For corrections like extinction, absorbtion and preffered orientation, it
might also want to know more about the experiment.

Least squares and other applications - like entropy maximisation and 
so on need to know the weight matrix coming from the data and the 
current calculated structure factors (and model weights).
Trying to make the various interfaces clean means that the least squares
should not know anything about crystallography.

The initial plan looks as if the program is going to trade memory
for complexity. 
That is to say, that very large intermediate results could be generated
in order to avoid different parts of the program from knowing too much
about each other.
Could be inefficient (memory traffic can get slow), but the clarity 
is expected to be enough of a benefit to make this worthwhile.
Current trends in computers (this laptop has 1GB of RAM) suggest
that running out of memory will not be a problem unless we do 
something really stupid.

\section{Essential variables}

For the peak fitting program we need a list of peaks. 
Each peak must define limits, height and a peak shape.
Optionally it might have other properties (like hkl, crystal phase),
but these are not interesting for the fitting program.
It is only going to compute the pattern and derivatives
with respect to the intensity, position and shape parameters.
Introducing a peak into a pattern requires a peak shape 
function, which can convert a peak into pixels.
Rather than list all of the peaks and give a function for each,
we can introduce peak shape functions and give each function a list
of peaks.
The peak shape function is then the link between the crystal structure
and the data. 
The crystal structure will define a peak shape function and that function
is given to the pattern. 

So the thing to pass is a ``function''.

That function is going to need arguments, which will need names that a human can 
interpret.
It will have parameters which it should compute derivatives for, and
parameters which it will not compute derivatives for. 
	The least squares thing, whatever it is,  needs to know about those parameters,
but only enough to know what is being refined.

So the generic object is a ``peak function''.
It needs a crystal structure to know where the peaks are.
It needs some instrumental parameters to compute itself.
It can put a peak into the dataset.
It can supply derivatives of each peak with respect to
intensity, position and shape parameters (however many there are).


\section{Can you give me a reference for that?}

Linking data structures and functions together without copying
is one of the essential things we need to be able to do in order
to manage to build a big complex program but without making a 
dirty mess.
The reason for choosing python is that the author figured out
how to use that language to accomplish that actual task.
It does not prevent you from messing things up, but it does have
the magical property of being some kind of a pointer manipulator
and calculator, without actually ever mentioning things like
allocate, free, dangling or null.
An example should show what I mean.

@o reference.py
@{
""" 
A short script which just shows something about python and references

Demonstrates something about classes, objects and how we plan to link
various things together without copying
"""
from Numeric import *
class myclass:
   """ Class to hold some data of some form """
   def __init__(self, data):
      """ Constructor - saves data in object """
      self.data=data

class myfunction:
   """ Class for functions - they need data to act on """
   def __init__(self, object_to_get_data_from):
      """ Remember that data """
      self.object_to_get_data_from=object_to_get_data_from
   def __call__(self):
      """ Calling the function object just returns the data for now """
      return self.object_to_get_data_from.data

data_object_10 = myclass(array((10,10)))# Instantiate an object of type myclass
data_object_a  = myclass(array((3,3)))  # Instantiate another one
func_10 = myfunction(data_object_10)    # Instantiate myfunction which needs an
                                        # object arg having a member ".data"
func_a = myfunction(data_object_a)      # One for a
print "\\begin{verbatim}"               # So we can put output into document
print "This should print [10 10]", func_10() # Call the func_10
print "This should print [3 3]", func_a()    # Similarly for func_a
print "Now we modify the data objects and call the functions again"
data_object_10.data = data_object_a.data
data_object_a.data = array((5,5))
print "This previously printed [10 10]", func_10() # Call the func_10
print "This previously printed [3 3]", func_a()   # Similarly for func_a
print "So everything is an object"
print "Assignment of x=y means that x now points at what y is pointing at"
print "\\end{verbatim}" # So we can put output into document
@}

The output from that program is:

\input{reference.py.out}

Based on this output, some ideas start to crystallise.
The crystal structure (objects) can exist as entities knowing
nothing about diffraction data.
They will know how to make lists of peaks, lists
of potential variables and derivatives of peaks with
respect to the variables.
The data histograms will know nothing about crystals, only that
they are composed of the sums of components.
The components will be given by a peak shape function, which links
itself to a crystal structure.

Therefore there will be as many crystal structure objects as we
want and as many histogram objects as we want. 
Each time a crystal structure is placed into a histogram we
do this using a peak shape function, which can look at the peak
list from the crystal and the data from the histogram.
It would generate a computed dataset and derivatives of the 
goodness of fit indicator.
Further pseudo-observations such as restraints are then 
just other kinds of ``peak shape function'' which look
at the crystal structure and some data of their own (bond lengths
and angles etc).


\section{Sparse matrix structures}

A lot of the computations involve very large arrays, which
have mostly zero entries.
We need to be able to pass sparse matrices around working
out which are the non-zero entries.
Following (for now) the matrix structure in PRODD,
we put the parameters into a particular order and 
define the length of the row needed for 
a particular entry and where it starts and ends.

For now I can only think of rectangular non-symmetric
matrices and square symmetric matrices that are likely
to be really useful. 
The rectangular non-symmetric ones are for passing 
derivatives of peak lists with respect to intensity
and position.



@o sparsesquarematrix.py
@{
""" Sparse matrix data structure, mimic mprodd """
@< pycopyright @>
from Numeric import *
class sparsesquarematrix:
   def __init__(self,pointers=None,elements=None):
      self.pointers=pointers
      self.elements=elements     
   def makefromarray(self,a):
      """ Fill out the pointers and elements from array a """
      pass
@}

@o sparserectmatrix.py
@{
""" Sparse rectangular matrix """
@< pycopyright @>
from Numeric import *
class sparserectmatrix:
   def __init__(self,ntotal,ndense):
       self.densrows=zeros(ndense,Float32)
       self.sparse=zeros(ntotal-nsparse,Float32)
       self.ntotal=ntotal
   def dotproduct(self,i,vector):
       """Compute mat[:,i] .dot. vec """
       # sparse part
@}


The output file format from prodd comes from the following code:
It is not very nice, but seems to work OK.
   
\begin{verbatim}
      WRITE(IOUT1,2000)ITOT,LVARB
2000  FORMAT(I8,' Intensities and ',I8,' least squares variables')
      WRITE(IOUT1,2001)
2001  FORMAT(' Pointers and variables for the matrix')
      DO I=1,ITOT
       WRITE(IOUT1,'(2(I8,1X),A3,1X,3(I5,1X),E15.5)')
     &              I,IPOINT(I),'HKL',NINT(REFH(1:3,I)),F4PAR(1,I)
      ENDDO
      DO I=ITOT+1,ITOT+LVARB
        WRITE(IOUT1,'(2(I8,1X),A3)')I,IPOINT(I),'LSQ'
      ENDDO
      WRITE(IOUT1,2002)
2002  FORMAT(' Now comes the matrix, one point on each line')
      DO J=1,IPOINT(ITOT+LVARB)
        WRITE(IOUT1,2003)ASS(J)
2003    FORMAT(E15.5)
      ENDDO
100   RETURN
      END
\end{verbatim}


\section{Refinement}

The interfaces which model and data objects need 
to present to be used in refinement are going
to depend on what the refinement needs to know.
We want to do this in a flexible high-level way, but
at the same time, we want to do it efficiently.
The model must supply a list of parameters
which are refineable, and it must be able
to give derivatives with respect to data.
The data just needs to give observed values
and weights, where the weight matrix might
be diagonal or sparse non diagonal.

For refinement against a conventional $\chi^2_R$
(pattern fitting) we expect derivatives with
respect to dense row parameters and also intensity-like
parameters which are defined on a per peak basis.

For refinement against a non-diagonal weight matrix
we expect derivatives which are dense rows (appear
for all observations, like atomic positions) and
also those which are sparse (reflections flagged for 
R$_{free}$). 

So the refinement algorithm generally expects refinable
parameters which contribute to all observations, or
only some observations.

What about using pseudo-data for making nicer difference 
fouriers? In this case we would be refining structural 
parameters (contribute to all peaks) and adding a contribution
to each peak representing the unlocated part of the structure.
That contribution would have a particular average size
depending on the amount of missing density and be refined to
make the misfit between model and data be larger than
given by least squares (larger by an amount expected based
on Wilson type contributions of the unlocated part).
The least squares problem is then like the pattern fitting
problem, with the intensities of the peaks being heavily
restrained to follow a particular function (depending on
the $F_{calc}$ values and amount of missing density.

So in general - we need the number of ``dense'' parameters
and the number of ``sparse parameters'' coming from a model.
We also need to link up observations to the model, 
to be able to ask what is the derivative of this 
observation with respect to that parameter.
Probably the model itself should take care of that.

So the least squares finally needs the obs and calc
values (differences), a weight matrix, and the derivative
of each difference with respect to the parameters in the model.
The parameters will be split into a ``dense'' set, which
supply derivatives to everything, and a ``sparse set'' which
supply derivatives to only some observations.
So an array of derivatives per observation for the dense
parameters. 
Also an array of derivatives of sparse parameters, which
also identify the observation(s) they contribute to.

Programmatically we need to add in contributions to 
the least squares matrix and rhs using 
the equations in the introduction.
Making the sparse matrix for the DILS refinements
puts the dense rows at the end.
Using this method outlined above, the model is 
computed and we need to note which peak intensities
contribute to each data point.
Then we would loop over observations working 
out the non-zero pattern for the matrix.
We will assume that the model sorts any sparse parameters
it contains into the most sensible order first.

\subsection{An example of optimisation}

The details of how the optimisation will work are not
important just yet. 
In writing this text, a short program to fit a polynomial
was created to get an idea of the things we 
need to do to define the interface. Here it is:

@O denseleastsquares.py
@{
""" Least squares optimiser object """
@< pycopyright @>
from Numeric import *
from LinearAlgebra import generalized_inverse, inverse

class denseleastsquares:
   def __init__(self,nv):
      """ nv is the number of dense variables """
      self.lsqA=zeros((nv,nv),Float)
      self.rhs =zeros((nv)   ,Float)
      self.chi =0.0
      self.nobs=0
      self.nv  =nv
   def addcontribution(self,obs,calc,wt,g):
      """ 
      obs = obs value
      calc= calc value
      wt  = statistical weight (probably 1/sigma**2)
      g   = gradients ( d(calc)/d x ) x in range(nv)
      """
      self.nobs+=1
      self.chi +=(obs-calc)*(obs-calc)*wt
      for i in range(self.nv):
         self.rhs[i]=self.rhs[i]+g[i]*(obs-calc)*wt
         for j in range(self.nv):
            self.lsqA[i,j]=self.lsqA[i,j]+g[i]*g[j]*wt
   def solve(self,damping=None):
      """
      Find a solution to the lsq problem please
      damping is for optional Marquardt style damping
       (lsqA[i,i]*=damping implies damping >= 1 )
      """
      if damping != None:
         print "Damping = ",damping
         for i in range(self.nv):
            self.lsqA[i,i]*=damping
      scales=self.normalise()
      self.invA=inverse(self.norm)
      self.shifts=matrixmultiply(self.invA,self.nr)
      self.shifts=self.shifts/scales
      # Check if they solve the problem
      leftover=matrixmultiply(self.lsqA,self.shifts)-self.rhs
      moreshift=matrixmultiply(self.invA,leftover/scales)/scales
      self.shifts=self.shifts+moreshift
      # no chi^2 correction to error bars here!
   def clean(self):
      """ Wipe clean ready for new problem """
      self.lsqA[:,:]=0.0
      self.rhs[:]=0.0
      self.nobs=0
      self.chi=0.0
      # leave invA and esds
   def normalise(self):
      s=zeros(self.nv,Float)
      self.norm=self.lsqA.copy()
      self.norm.savespace(1)
      self.nr=self.rhs.copy()
      for i in range(self.nv):
         s[i]=sqrt(self.lsqA[i,i])
         self.norm[:,i]=self.norm[:,i]/s[i]
         self.norm[i,:]=self.norm[i,:]/s[i]
         self.nr[i]=self.nr[i]/s[i] 
      return s

@}

In order to give a real concrete and useful example
we want to fit a real model to a real dataset.
In this case a degree n polynomial will do for the 
model. 
The data can be just some random numbers.
Eventually this could make a handy background function, 
but for the numerical properties of very high order
polynomials.

@o polymodel.py
@{
""" Computes a polynomial model to x,y data """
@< pycopyright @>
from Numeric import *
class polymodel:
   def __init__(self,degree,coefficients):
      """ Degree = number of coefficients + 1 
          since zeroth order is a constant """
      self.degree=degree
      if coefficients.shape != (degree+1,):
         print coefficients.shape,(degree+1,)
         raise Exception("Wrong number of coefficients")
      self.ai=coefficients
   def compute(self,x):
      """ x - dependent variable val , compute gradients too """
      y=self.ai[0]
      g=ones(self.degree+1,Float)
      for i in range(1,self.degree+1):
         y   = y + self.ai[i]*pow(x,i) 
         g[i]= pow(x,i)
      return y,g

if __name__=="__main__":
   p=polymodel(3,array([.1,.2,.3,.4],Float32))
   y,g = p.compute(3)
   print y,g
@}
 
Now make a little program which generates the data and
fits the polymodel.

@O testdenseleastquares.py
@{
""" Fit a function to some data """
@< pycopyright @>

from denseleastsquares import *
from polymodel import *
from Numeric import *
from RandomArray import random
# Make x data
np=100
nv=20
x=(arange(np)-np/2)/(np*1.)
start=random(nv)
for i in range(nv): start[i] *= (nv-i)
print "Start = ",start
p=polymodel(nv-1,start)

y=zeros(np,Float)

for i in xrange(y.shape[0]):
   y[i]=p.compute(x[i])[0]

# Data are now in x,y
model=zeros(nv,Float)
q=polymodel(nv-1,model)

lsq=denseleastsquares(nv)

for j in range(20):
   for i in xrange(y.shape[0]):
      yc,g=q.compute(x[i])
      lsq.addcontribution(y[i],yc,y[i],g)

   lsq.solve()
   es=sqrt(sum((model-start)*(model-start)))
   for i in range(nv): model[i]=model[i]+lsq.shifts[i]
   print "Chi^2=%g, nobs=%d, Reduced Chi^2=%g"%(lsq.chi,lsq.nobs,lsq.chi/lsq.nobs),
   print "error sum = %g" % (es)
   lsq.clean()
@}


Testing out that script debugged the denseleastsquares class.
It improved a huge amount when double precision replaced
single precision. 
Whether it makes much difference if the matrix is normalised
before inversion is not so clear.
Things to do now are to put the singular value decomposition
in to find out how that works (a later chapter on optimisation
methinks).

Strange behaviour with the polynomial function - it becomes very
badly behaved once the order becomes higher. 
Perhaps it is worth investigating orthogonal polynomials
and their properties.
This would be crappy as a background function as it is possible
to get a very good $\chi^2$ with a very poor reproduction
of the original parameters.


\subsection{An interface to optimisation}

The polynomial fitting example gives an example of 
a full program in python which fits a mathematical model
to some data.
We can have a look at a trimmed down version of 
the optimisation object created:
\begin{verbatim}
class denseleastsquares:
   def __init__(self,nv):
      """ nv is the number of dense variables """
   def addcontribution(self,obs,calc,wt,g):
      """ 
      obs = obs value
      calc= calc value
      wt  = statistical weight (probably 1/sigma**2)
      g   = gradients ( d(calc)/d x ) x in range(nv)
      """
   def solve(self,damping=None):
      self.shifts
   def clean(self):
   def normalise(self):
\end{verbatim}
It is a python class, so we could have several instances of the 
optimiser present in the program at the same time, for comparing
different models without losing any information.
Before it could get started it needed to know how
many variables were being refined (\code{nv}).
Each observation is added into the problem using the \code{addcontribution}
function and set of parameter shifts is found from the \code{solve}
function.
It also offers a \code{clean} function for removing all the previously added
contributions and getting ready for another refinement cycle and also
a normalise function which is only used internally.

While this might be OK for fitting a polynomial we need something more
generalised for data where the weight is a non-diagonal matrix and 
the gradients are frequently a sparse vector.

Dealing with a sparse vector could be done via a sparse vector object.
This will have a number of elements and arrays of pointers which say 
which vector element is which and what the value is.
We also need to find a generalised way of dealing with constraints,
so in some senses it makes sense to have a potentially large array 
of derivatives with respect to parameters which are converted
into an array of derivatives with respect to least squares variables.
From the point of view of the optimiser, this is playing around 
with the model and has nothing to do with optimising things, so we
should deal with the parameter to variable transformations outside
of this class in building our refinement model.

Keeping the least squares simple means leaving the model building out.
Deciding how to propagate errors into and out of the optimisation 
routines will be done elsewhere. 
This allows different optimisers to be used in ways which are easier to 
deal with,
The optimiser above assumes there is no covariance between observations.
When there is covariance we need derivatives of the correlated parameters
at the same time.
A simple approach is to make all of the derivatives into a matrix
and compute the matrix products ($[\alpha] = \mb{MAM}$.
This is heavy on memory but we can use the sparse structure in $\mb{A}$ 
for saving on computations. 
So in this way we would have to make a matrix of gradients which are 
full of $\partial y_obs/\partial x$.
For pattern fitting that makes  a potentially very large matrix, which 
should only be storing the elements which are needed (derivatives of
peaks are zero in a lot of places).

In practice we have sparse derivative vector for pattern fitting 
and a sparse weight matrix for structure refinement.
We might deal with these two case separately, rather than trying
to get something working immediately with sparse vectors and 
matrices at the same time. 
The sparse weight matrix perhaps needs to offer it's
own optimisation scheme.
With a set of sparse derivative vectors we need to determine
the matrix structure which will be used.




\section{Crystal structures}

The minimal thing we need is a python object which can
supply the information needed about a particular
crystal structure.
Information which will be needed is the following:
\begin{itemize}
\item List of hkl peaks
\item Computed peak intensities
\item Variable parameters
\item Derivatives of peaks with respect to variables
\item Optionally a list of pseudo-observations (restraints)
\end{itemize}
The crystal structure may or may not contain atoms - it
might represent some kind of a density (for example, modelling
a C$_60$ molecule by a sphere of scattering density).
The simplest kind of crystal structure, which all other
crystal structures must be able to mimic (inherit from), is
a list of peaks.
This kind of ``crystal structure'' would only be used for 
peak fitting prior to indexing.

@o peak_array.py
@{
""" 
Base class for crystal structures
Defines interface to the rest of linarp
"""
@< pycopyright @>
from Numeric import *
class peak_array:
   """
   Minimally holds a Numeric array of peak positions (in Q=1/d**2)
   and intensities (in arbitrary units)
   """
   def __init__(self,positions=None,intensities=None):
      """ Absolute minimum is peak list in Q """
      self.positions=positions
      self.intensities=intensities
   def test(self):
      if self.positions != None and self.intensities != None:
         if self.positions.shape != self.intensities.shape:
            raise Exception("Corrupted peak array")   
      if self.positions != None and self.intensities == None:
            raise Exception("Corrupted peak array")   
      if self.positions == None and self.intensities != None:
            raise Exception("Corrupted peak array")   
@}


For peak fitting we need each peak to have a shape, intensity
and width. 
The width will come from the peak shape function, which links
a crystal structure to a model.
The intensity will come from the crystal structure. 
It might initially be zero, and the crystal structure will
either let us update it when peak fitting by supplying
it directly as a variable, or it will update the value
itself by applying shifts to atomic co-ordinates.
For crystal structure refinement we need the peaks to also have
hkl indices, intensities and derivatives of these with
respect to any parameters.


\subsection{Peak fitting type crystal structures}

These are a kind of crystal which is used for indexing.
It doesn't have a unit cell or spacegroup, and only
records peak positions and intensities.
To use it for peak fitting it will need to provide
derivatives of positions and intensities and
the number of refinable parameters.

We immediately find out that we need to work out
how to deal with sparse matrices in a clever way here.
The peak fitting object is going to want to provide
derivatives with respect to position and intensity
for each peak, but the derivatives of the other 
peaks are clearly zero.
If a unit cell is defined then all of the peaks will
have derivatives on position with respect to 
unit cell parameters.
Therefore the derivatives could be partitioned into
a sparse and dense set.
Something for passing sparse matrices around is needed.



@o pf_xtal.py
@{
""" 
Dummy structure for peak fitting
"""
@< pycopyright @>
from Numeric import *
class pf_xtal(peak_array):
   def addpeaks(self,positions,intensities):
      """ Append further peaks """
      if self.positions==None:
         self.positions=positions
         self.intensities=intensities
      else:
         newsize=positions.shape[0]+self.positions.shape
         self.positions=resize(self.positions,newsize)
         self.positions[-positions.shape[0]:]=positions
   def n_i_pars(self):
      """ Return the number of refinable intensity parameters 
      """
      return self.intensities.shape[0]
   def derivatives(self):
      """ Return the derivative of peaks w.r.t the
          refinable parameters """
      
@}




\chapter{Input data formats}

Reading raw datafiles from the many varied instruments which exist
is a cause for serious pain for anyone who has done a lot 
of powder diffraction.
We'll make a data object here which we hope will be general
enough for whatever formats come along.
For most people an array of ``x'' values combined with arrays
of intensity and error bar are enough. 
For some people, especially those using area detectors, we 
might want to let the ``x'' array have more than one number 
in it (position on detector and maybe sample setting angle).

Other items which are needed about the experimental setup
can all be stored in a dictionary. For example, wavelengths,
monochromator angles, polarisation, temperature, pressure, phase
of the moon, etc. 

The generic data object holds an x array, a y array and and e array 
and a dataitems dictionary. For a 2D image you can make the arrays 1D
still, but put a stride in the dictionary. 
Depends on personal taste I guess?

@o powderdata.py
@{
""" Dataset object for powder data """
@<pycopyright@>
from Numeric import *
class powderdata:
   def __init__(self,x,y,e,d):
      """ x = x values, y = y values, e = errors, d = dictionary """
      self.x=self.xfull=x
      self.y=self.yfull=y
      self.e=self.efull=e
      self.w=self.wfull=(1./(self.efull*self.efull)).astype(Float32)
      self.d=d
      self.npoints=self.x.shape[0]

   def setrange(self,lh=None):
      """ Set the active range """
      if lh==None:
         self.x=self.xfull
         self.y=self.yfull
         self.e=self.efull
         self.w=self.wfull
         self.d["active_range"]=None
         self.npoints=self.x.shape[0]
      else:
         l=searchsorted(self.xfull,lh[0])
         h=searchsorted(self.xfull,lh[1])
         self.x=self.xfull[l:h]
         self.y=self.yfull[l:h]
         self.e=self.efull[l:h]
         self.w=self.wfull[l:h]
         self.d["active_range"]=lh
         self.npoints=self.x.shape[0]
  
   def weightmatvec(self,vector):
      """ 
      Converts d{yalc}/d{parameter} into d{chi^2}/d{parameter} via 1/sigma^2
      """
      return vector*self.w
@}

The various dictionary items which you might need are defined here.
Please add to this list in the future or make your file reader
use the things in here already.
Other parts of the code are going to rely on finding the information
they need in this dictionary.

\begin{itemize}
\item wavelength - for constant wavelength data
\item zero -  if there is a zero shift to add to x
\item pathlength  - for ToF neutron data
\item angle - for ToF neutron of energy dispersive x-ray
\item polarisation - in range 0 to 1.
\item temperature
\item pressure
\item time
\item title
\item formula
\item radiation - x-ray, neutron, electron
\item xunits - these determine CW versus ToF versus ED
 \begin{itemize}
  \item degrees - implies 2$\theta$, constant wavelength
  \item keV - implies photon energy for energy dispersive
  \item microseconds - implies neutron ToF
 \end{itemize}
\item yunits
\item mur  - absorbtion coefficient (for cylinders)
\item mut - absorption coefficient (for plates)
\item comment - anything you like
\item fromfile - where the data were read in from
\item Some system for sample orientation angles - define please - FIXME
\end{itemize}


Finally - a handy little function for making filename from numbers
and stems and dealing with multiple datasets:

@d makename
@{
def makename(stem,num,extn,width=4):
   """
   Generates filenames "stem"+nnnn+"extn" 
   eg: Bruker: sample101.0000 is "sample101."+0+"",width=4
       edf   : sample0000.edf is "sample"+0+".edf",width=4
       chi   : sample0000.chi is "sample"+0+".chi",width=4
       inp   : 0.inp is ""+0+".inp",width=0
   """
   if width>0:
      fs="%s%0"+"%d"%(width)+"d%s"
   else:
      fs="%s%d%s"
   return fs%(stem,num,extn)
@}

\section{Reading epf files}

At the esrf the usual dataformat for powders is just a list of 
x y error values.
We'll make a little object for reading the files.

@o epffile.py
@{
""" Epf file reading object """
@< pycopyright @>
from powderdata import powderdata
from Numeric import *
class epffile(powderdata):
   def __init__(self,filename,**kwargs):
      """ Reads a filename and hopes to get wavelength etc as optional args """
      x=[] ; y=[]; e=[]; d={}
      d["fromfile"]=filename   # defaults
      for line in open(filename).readlines():
         xye=map(float,line.split())
         x.append(xye[0])
         y.append(xye[1])
         e.append(xye[2])
      d["xunits"]="2theta"
      d["yunits"]="arbitrary"
      for key,val in kwargs.items(): d[key]=val
      powderdata.__init__(self,array(x),array(y),array(e),d)
@}

\section{Reading GSAS files}

Something should go here for reading GSAS files

\section{Reading powbase files}

Something to read the powbase datafile format

From http://www.cristal.org/powbase/add.html
\begin{itemize}
\item Comments beginning with ``\#'' or ``!''
\item Three numbers which $2\theta$ start, step and end.
\item Raw intensities in free format
\end{itemize}

@o powbase.py
@{
""" Powbase file reading object """
@< pycopyright @>
from powderdata import powderdata
from Numeric import *
class powbase(powderdata):
   def __init__(self,filename,**kwargs):
      y=[]; e=[]; d={}
      d["fromfile"]=filename
      low=None
      for line in open(filename).readlines():
         if line[0] in ["!","#"]:
            if d.haskey("comment"):
               d["comment"]="%s\n%s"%(d["comment"],line)
            else:
               d["comment"]=line
         elif low==None:
            # read tth low, step, high
            [low,step,high]=map(float,line.split())
            x=arange(low,high+step,step)
         else:
            vals=map(float,line.split())
            for v in vals: y.append(v)
      t=array(y,Float)
      e=sqrt(y)
      d["wavelength"]=1.54036
      d["xunits"]="2theta"
      d["yunits"]="Counts"
      for key,val in kwargs.items(): d[key]=val
      powderdata.__init__(self,array(x),array(y),array(e),d)
@}

\section{Reading MCA files}

This is an ESRF spec format with an MCA data array.

@o mcadata.py
@{
""" MCA esrf spec file reader """
@< pycopyright @>
from powderdata import powderdata
from Numeric import *
class mcadata(powderdata):
   def __init__(self,filename,**kwargs):
       y=[] ; e=[] ; d={}
       d['fromfile']=filename 
       on=0
       for line in open(filename).readlines():
          if on==1:
             nums=map(int,(line.split("\\")[0]).split())
             y.extend(nums)
          if line[0]=="#":
             d[line.split()[0]]=line
          if line[0]=="@@":
             on=1
             l=line[2:].split("\\")[0]
             nums=map(int,l.split() )
             y.extend(nums)
       x=array(range(len(y)))
       y=array(y)
       e=sqrt(y+0.5) # Bayesian 0.5 to avoid zero errorbars
       powderdata.__init__(self,x,y,e,d)

 
if __name__=="__main__":
   import sys
   obj=mcadata(sys.argv[1])
   

@}

\section{Reading cif (data) files}

Read in the powder data related parts of cif files.
Other things (like crystal structures) to go elsewhere.

\section{Reading 2 dimensional images}

There are many possible image formats which we might eventually
like to be able to read. For now something is included for a 
Bruker Smart area detector and the ESRF EDF data file format.
Neither of these decoders is a full implementation, but intended
to be just enough for a quick fix for reading data for beamline
ID11 at the ESRF.

For most image reading applications we are supplied with an
array dumped as a stream of bytes. Conversion into a python
array of that stream of bytes is given by the function below.

@d readbytestream
@{
def readbytestream(file,offset,x,y,nbytespp,datatype='int',signed='n',
                   swap='n',typeout=UInt16):
   """
   Reads in a bytestream from a file (which may be a string indicating
   a filename, or an already opened file (should be "rb"))
   offset is the position (in bytes) where the pixel data start
   nbytespp = number of bytes per pixel
   type can be int or float (4 bytes pp) or double (8 bytes pp)
   signed: normally signed data 'y', but 'n' to try to get back the right numbers
     when unsigned data are converted to signed (python has no unsigned numeric types.)
   swap, normally do not bother, but 'y' to swap bytes
   typeout is the Numeric type to output, normally UInt16, but more if overflows occurred
   x and y are the pixel dimensions
   """
   tin="dunno"
   len=nbytespp*x*y # bytes per pixel times number of pixels
   if datatype=='int' and signed=='n':
      if nbytespp==1 : tin=UInt8
      if nbytespp==2 : tin=UInt16
      if nbytespp==4 : tin=UInt32
   if datatype=='int' and signed=='y':
      if nbytespp==1 : tin=Int8
      if nbytespp==2 : tin=Int16
      if nbytespp==4 : tin=Int32
   if datatype=='float':
      tin=Float32
   if datatype=='double' :
      tin=Float64
   if tin=="dunno" :
      raise SyntaxError, "Did not understand what type to try to read"
   opened=0
   if type(tin) == type(file):  # Did we get a string or a file pointer?
      f=open(file,'rb')
      opened=1
   else:
      f=file
   f.seek(offset)
   if swap=='y':
      ar=array(reshape(byteswapped(fromstring(f.read(len),tin)),(x,y)),typeout)
   else:   
      ar=array(reshape(fromstring(f.read(len),tin),(x,y)),typeout)
   if(opened):f.close()
   return(ar)
@}


\subsection{Bruker area detector format}

We begin with something to read in the header section of the file:

@d readbrukerheader
@{
def readbrukerheader(file):
   """
   Reads a Bruker file header into a Python dictionary
   file=filename or file pointer
   """
   s="string"                 # s is a string
   if type(s) == type(file):  # if arg is a string, open file, else treat as file object
      f=open(file,"rb")
      opened=1
   else:
      f=file
      opened=0                # opened var flags to close again if we open the file
   i=80
   block=f.read(512)          # always start with a 512 byte header
   Header={}                  # dict to take results
   while i < 512 :            # wander along the 512 bytes
      key,val=block[i-80:i].split(":",1)   # uses 80 char lines in key : value format
      key=key.strip()         # remove the whitespace (why?)
      val=val.strip()
      if Header.has_key(key):             # append lines if key already there
         Header[key]=Header[key]+'\n'+val
      else:
         Header[key]=val
      i=i+80                  # next 80 characters
   nhdrblks=int(Header['HDRBLKS'])    # we must have read this in the first 512 bytes.
   # print "got first chunk, headerblocks is",nhdrblks
   # Now read in the rest of the header blocks, appending to what we have
   block = block[i-80:512] + f.read(512*(nhdrblks-1))
   j=512*nhdrblks
   while i < j :
      # print i,"*",block[i-80:i].strip(),"*"
      if block[i-80:i].find(":") > 0:          # as for first 512 bytes of header
         key,val=block[i-80:i].split(":",1)    
         key=key.strip()
         val=val.strip()
         if Header.has_key(key):
            Header[key]=Header[key]+'\n'+val
         else:
            Header[key]=val
      i=i+80
   Header['datastart']=f.tell()                # make a header item called "datastart"
   if(opened):f.close()
   return Header     # s
@}

And now the bulk of the data - doing something about the overlaps at the end
of the file.

@d readbruker
@{
@< readbrukerheader @>

def readbruker(file):
   """
   Reads in a Bruker file, returning the data and header
   file may be a string or file object
   FIXME we should later modify to take ROI ranges somehow (xmin,xmax,ymin,ymax)
   """
   s="string"                 # s is a string
   if type(s) == type(file):  # if arg is a string, open file, else treat as file object
      f=open(file,"rb")
      opened=1
   else:
      f=file
      opened=0       
   Header=readbrukerheader(f) # pass in the file pointer, it stays open
   npixelb=int(Header['NPIXELB'])   # you had to read the Bruker docs to know this!
   rows   =int(Header['NROWS'])
   cols   =int(Header['NCOLS'])
   # We are now at the start of the image - assuming readbrukerheader worked
   size=rows*cols*npixelb
   data=readbytestream(f,f.tell(),rows,cols,npixelb,datatype="int",signed='n',swap='n')
   no=int(Header['NOVERFL'])        # now process the overflows
   if no>0:   # Read in the overflows
       # need at least Int32 sized data I guess - can reach 2^21
       data=data.astype(UInt32)
       # 16 character overflows, 9 characters of intensity, 7 character position
       for i in range(no):
          ov=f.read(16)
          intensity=int(ov[0:9])
          position=int(ov[9:16])
          r=position%rows           # relies on python style modulo being always +
          c=position/rows           # relies on truncation down
          #print "Overflow ",r,c,intensity,position,data[r,c],data[c,r]
          data[c,r]=intensity
   f.close()
   Header["rows"]=rows
   Header["columns"]=cols
   return Header,data
@}


A little script to test out our bruker reading capabilities...

@o testbruker.py
@{
from Numeric import *

import imagereaders

if __name__ == "__main__":
   import sys
   if len(sys.argv) < 2:
      print "Usage: ",sys.argv[0]," filename"
      sys.exit()
   file=sys.argv[1]
   h,d=readbruker(file)  # the command to get the file read in
   for k,v in h.items(): # print out the header info
     print k,":",v
   print "Data shape:",d.shape # image shape
   print "Data max, min, sum:",maximum.reduce(ravel(d)),   \
                              minimum.reduce(ravel(d)),   \
                              add.reduce(ravel(d))
@}


This covers the basics of reading Bruker files. They are going to come
back as Unsigned Integer 8, 16 or 32 bit types. Optimisations to only
read in certain parts of the image for doing faster computations
of rocking curves might be nice, but that doesn't seem so important.
The work involved in doing it seems a waste of time, as if you are
going to be there looking for a peak to do a rocking curve of, you've
already read the entire image in to choose it. Only really worthwhile
for 8mb images. FIXME (do this getroi crap sometime for readbytestream)


\subsection{Read ESRF edf files}

We implement a very brief edf file reader here. The format is actually
more extensive - you can place multiple images in one file etc. \emph{However}
for reasons which are best not discussed this reader is only expected to
work for images created, incorrectly, at beamline ID11. There is a spare
pixel, which means you have to read back from the end of the files. Please
don't ask!

@d readid11edf
@{
def readid11edf(filename):
   f=open(filename,"rb")
   # Header is 1024 byte blocks, terminated by "}"
   fh=f.read(1024)
   i=1023
   while fh[-2:]!="}\n": 
      fh+=f.read(1024)
   # Interpret header
   headeritems=fh[1:-1].split(";")
   hd={}
   for item in headeritems: 
      if item.find("=")>0:
         hd[item.split("=")[0].lstrip().rstrip()]=item.split("=")[1].lstrip().rstrip()
   hd["rows"]=int(hd["Dim_1"])
   hd["columns"]=int(hd["Dim_2"])
   f.seek(-int(hd["Size"]),2) # go to position offset from end
                              # assumes last byte is last pixel (first is occasionally not!)
   datastring=f.read(int(hd["Size"]))
   f.close()
   # Convert datastring to Numeric array
   if hd["DataType"]=="UnsignedShort": numerictype=UInt16
   else: 
      raise TypeError("Unimplemented edf filetype")
   if hd["ByteOrder"]=="LowByteFirst": 
      ar=reshape(
            fromstring(datastring,numerictype),
            (hd["rows"],hd["columns"]) )
   else:
      ar=reshape(
            fromstring(datastring,numerictype).byteswapped(),
            (hd["rows"],hd["columns"]) )
   return hd,ar
@}

\subsection{Reading image files}

Putting together the bruker and edf readers we have:

@o imagereaders.py
@{
""" 
Image reading functions. They return a numeric array representing
the image and a dictionary of thing which were found in headers
"""
    
@< pycopyright @>

from Numeric import *

@< makename @>
@< readbytestream @>
@< readbruker @>
@< readid11edf @>
@}



\section{The correlated integrated intensity stuff from prodd}


A python object which can read in the sparse matrix structure
coming from prodd and do something useful with it. Eventually we need
to support the same kinds of interfaces as the powderdata object
for things like plotting but probably some extra stuff like
finding out what the hkls of the peaks are.


\section{The Fortran Library}

The fortran version of Cholesky decomposition is as follows. This
fortran was originally implmented in the prodd computer program
for carrynig out large scale intensity extractions.

@D fchol
@{
      SUBROUTINE CHOLDC(A,IP,N)
C
CDEC$ attributes alias : '_choldc_' :: CHOLDC     
C $
C Cholesky decomposition routine. Part of MATDIL above
C 
C No longer safe - it relies on the matrix being preconditioned to 
C give diagonal dominance (ie multiply all diagonal elements by 0.1)
C
C A IS THE MATRIX TO BE DECOMPOSED
C IP IS A SET OF POINTERS TO THE DIAGONAL ELEMENTS OF A
C N IS THE DIMENSION OF THE MATRIX A
C
C The decomposition overwrites the input matrix
C
C Essentially the same as the envelope method in George and Liu, "Computer
C Solution of Sparse Positive Definite Systems", Prentice Hall, 1981
C
C Note that the array indexing is different from that text
C
      REAL A(*)
      DOUBLE PRECISION SUM
      INTEGER N,IP(N)
C During this loop over 'I' we build up the Cholesky factor,
C Get L(1,1) as a special case
      A(1)=SQRT(A(1))
      DO I=2,N
C Solve the previous lower triangle system... only goes back to start of row
        IR=IP(I)-IP(I-1)-1
        IPI=IP(I)-I
C J=1 is also a special case (J-1 and I-1 don't exist)
        IF(I-IR.EQ.1) THEN
          A(IPI+1)=A(IPI+1)/A(1)
          JMIN=2
        ELSE
          JMIN=I-IR
        ENDIF
C Now solving for the row up to the element before diagonal element
        DO J=JMIN,I-1
          SUM=A(IPI+J)
          JR=IP(J)-IP(J-1)-1
          KM=MAX(J-JR,I-IR)
          IPJ=IP(J)-J
C This is the loop which would cost N^3 if not sparse
          DO K=J-1,KM,-1
            SUM=SUM-A(IPI+K)*A(IPJ+K)
          ENDDO
          A(IPI+J)=SUM/A(IP(J))
        ENDDO
C Compute Lii
        SUM=A(IP(I))
        DO K=I-1,I-IR,-1
          SUM=SUM-A(IPI+K)*A(IPI+K)
        ENDDO
C Guard against negative definite systems - sets sqrt(neg) to 
C be a large number - effectively adding a strong restraint to
C the thing that blew up.
7       IF(SUM.LT.0.0) THEN
          SUM=1.e10
        ENDIF
        A(IP(I))=SQRT(SUM)
      ENDDO
      RETURN
      END
@}


@d matvec
@{      subroutine MATVEC(A,IP,X,B,N)
C
CDEC$ attributes alias : '_matvec_' :: MATVEC
C $
C A is the matrix
C IP are the pointers to diagonal elements
C X is the vector
C B is set to equal A.x
C N is the dimension of the matrix
      integer N, IP(N)
      real X(N),B(N),A(*)
      DO I=1,N
       B(I)=0.
      ENDDO
      B(1)=B(1)+A(1)*X(1) ! first element
      DO I=2,N     
C Diagonal (to avoid counting it twice)
        B(I)=B(I)+X(I)*A(IP(I))
        DO J=IP(I-1)+1,IP(I)-1
          B(I)=B(I)+A(J)*X(I-IP(I)+J) ! Rows 
          B(I-IP(I)+J)=B(I-IP(I)+J)+A(J)*X(I) !  Columns
        enddo
      enddo
      end subroutine matvec
@}

And solving in fortran is done via:

@d solvec
@{
      SUBROUTINE SOLVE(A,IP,B,X,N)
C
CDEC$ attributes alias : '_solve_' :: SOLVE
C $ 
C
C A  = Cholesky decomposed matrix
C IP = Pointers to the diagonal elements
C B  = RHS of LSQ problem
C X  = Vector to contain solution vector
C N  = Dimension of the problem
      INTEGER IP(N),N
      REAL A(*),B(N),X(N)
C Matrix is decomposed, now solve L.Y=B and L^T.X=Y eg L (L^T.X) = B
      X(1)=B(1)/A(1)
      DO I=2,N
        SUM=B(I)
        IR=IP(I)-IP(I-1)-1
        DO K=I-1,I-IR,-1
C K is going back along a row
          SUM=SUM-A(IP(I)-I+K)*X(K)
        ENDDO 
5       X(I)=SUM/A(IP(I))
      ENDDO
      DO I=N,1,-1
        SUM=X(I)
        DO 6 K=I+1,N
C Now K is going vertically up the rows
C This loop could be better optimised....      
          IF(IP(K)-K+I.LE.IP(K-1))GOTO 6
          SUM=SUM-A(IP(K)-K+I)*X(K)
6       CONTINUE
        X(I)=SUM/A(IP(I))
      ENDDO
      RETURN
      END
@}

@o libchol.f
@{
@< fchol  @>
@< matvec @>
@< solvec @>
@}



@O pylibchol.c
@{
/* Python-Numeric wrapper for fortran cholesky routines */
#include <Python.h>
#include "Numeric/arrayobject.h"
typedef int* INTEGER;
typedef float*  REAL;

static char getrefcount_docstring[] =
"Returns the current reference count for object";

static PyObject * getrefcount (PyObject *self, PyObject *args){
PyObject *in=NULL;
if(!PyArg_ParseTuple(args,"O",&in)) return NULL;
return Py_BuildValue("i",in->ob_refcnt-1); /* the minus 1 is for this reference here */
}

   

static char fcholdc_docstring[] =
"L=fcholdc(A,IP,N)\n\n"\
"Compute the sparse Cholesky decomposition of a matrix, A.\n"\
"Pointers to the diagonal elements are stored in IP.\n"\
"N indicates the number of rows/columns of the matrix to consider (0:N)\n"\
"The result, L, is stored in the same way as A and indexed by IP[0:N]\n\n"\
"Assumes the problem is non-singular, with bizarre damping if it is not.\n\n"\
"A  : Sparsely stored matrix values\n"\
"IP : Pointers to diagonal elements (fortran indexing)\n"\
"N  : Integer - dimension of full square matrix A and vector N in Ax=b";

static PyObject * fcholdc (PyObject *self, PyObject *args){

   extern void choldc_( REAL, INTEGER, INTEGER);
   PyObject *array, *pointers;
   int dim[1],j,dimensions[1];
   PyArrayObject *matrix, *ip, *result;
   
   /* Read the arguments */
   if(!PyArg_ParseTuple(args,"OOi",&array,&pointers,dim))return NULL;
 
   matrix=(PyArrayObject *)
      PyArray_ContiguousFromObject(array,PyArray_FLOAT,1,1);

   if (matrix==NULL){
    PyErr_SetString(PyExc_ValueError, "Arg 1 did not make it to Float32 array");
    return NULL;
   }

   ip=(PyArrayObject *)
      PyArray_ContiguousFromObject(pointers,PyArray_INT,1,1);

   if (ip==NULL){
    PyErr_SetString(PyExc_ValueError, "Arg 2 did not make it to Int32 array");
    return NULL;
   }

   /* Check arguments, dim(IP)>=N and dim[A]>=IP[n] */

   if (ip->dimensions[0] < dim[0]) {
    PyErr_SetString(PyExc_ValueError, "Must have Arg 3 <= Arg2.shape[0]");
    return NULL;
   }
   j= (* (int *) (ip->data + (dim[0]-1) * ip->strides[0]));

   if (matrix->dimensions[0] < j) {
    PyErr_SetString(PyExc_ValueError, "Arg1 not big enough for arg2[arg3]");
    return NULL;
   }
   dimensions[0]=j;
   result=(PyArrayObject *)PyArray_FromDims(1,dimensions,PyArray_FLOAT);
   if (result==NULL){
    PyErr_SetString(PyExc_ValueError, "Wow, could not make a results array");
   }
   /* choldc overwrites in place, so we copy in matrix */
   for(j=0;j<dimensions[0];j++){
      (* (float *) (result->data + j * result->strides[0])) =
      (* (float *) (matrix->data + j * matrix->strides[0])); }
   choldc_((REAL) result->data, (INTEGER) ip->data, dim);

   PyArray_XDECREF(matrix);
   Py_XDECREF(array);
   PyArray_XDECREF(ip);
   PyArray_XDECREF(result);
   return PyArray_Return(result);
}

static char fmatvec_docstring [] =
"b=fmatvec(A,IP,X,N)\n"\
"\nComputes the matrix vector product Ax with A in sparse storage\n"\
"indexed by IP and using only up to dimension N of vector and\n"\
"matrix.\n\n"\
"A  : Sparsely stored matrix values\n"\
"IP : Pointers to diagonal elements (fortran indexing)\n"\
"X  : 1-D vector\n"\
"N  : Integer - dimension of full square matrix A and vector N in Ax=b";

static PyObject * fmatvec (PyObject *self, PyObject *args){

   extern void matvec_( REAL, INTEGER, REAL, REAL, INTEGER);
   int dim[1],j;
   PyArrayObject *matrix, *ip, *x, *result;
   
   /* Read the arguments */
   if(!PyArg_ParseTuple(args,"O!O!O!i",
                       &PyArray_Type, &matrix,
                       &PyArray_Type, &ip,
                       &PyArray_Type, &x ,dim))return NULL;
   /* Check arguments, dim(IP)>=N and dim[A]>=IP[n] */

   if (ip->dimensions[0] < dim[0]) {
    PyErr_SetString(PyExc_ValueError, "Must have Arg 4 <= Arg2.shape[0]");
    return NULL;
   }
   if (x->dimensions[0] < dim[0]) {
    PyErr_SetString(PyExc_ValueError, "Must have Arg 4 <= Arg3.shape[0]");
    return NULL;
   }
   j= (* (int *) (ip->data + (dim[0]-1) * ip->strides[0]));

   if (matrix->dimensions[0] < j) {
    PyErr_SetString(PyExc_ValueError, "Arg1 not big enough for arg2[arg4]");
    return NULL;
   }
   /* Check types : expecting 4 byte reals and integers */


   /* Make results */
   result=(PyArrayObject *)PyArray_FromDims(1,dim,PyArray_FLOAT);
   if (result==NULL){
    PyErr_SetString(PyExc_ValueError, "Wow, could not make a results array");
   }

   matvec_((REAL) matrix->data,
           (INTEGER) ip->data, 
           (REAL) x->data, 
           (REAL) result->data, 
            dim);
   PyArray_XDECREF(matrix);
   PyArray_XDECREF(ip);
   PyArray_XDECREF(result);
   PyArray_XDECREF(x);
   return PyArray_Return(result);
}





static char fsolve_docstring [] =
"b=fsolve(A,IP,B,N)\n"\
"\nCSolves the system (L L^T)x=b with L in sparse storage n"\
"indexed by IP and using only up to dimension N of vector and\n"\
"matrix.\n\n"\
"L  : Sparsely stored matrix values (from fcholdc)\n"\
"IP : Pointers to diagonal elements (fortran indexing)\n"\
"B  : 1-D vector\n"\
"N  : Integer - dimension of full square matrix A and vector N in Ax=b";

static PyObject * fsolve (PyObject *self, PyObject *args){

   extern void  solve_( REAL, INTEGER, REAL, REAL, INTEGER);
   PyObject *array, *pointers, *vec;
   int dim[1],j;
   PyArrayObject *matrix, *ip, *x, *result;
   
   /* Read the arguments */
   if(!PyArg_ParseTuple(args,"OOOi",&array,&pointers,&vec,dim))return NULL;
 
   matrix=(PyArrayObject *)
      PyArray_ContiguousFromObject(array,PyArray_FLOAT,1,1);

   if (matrix==NULL){
    PyErr_SetString(PyExc_ValueError, "Arg 1 did not make it to Float32 array");
    return NULL;
   }

   ip=(PyArrayObject *)
      PyArray_ContiguousFromObject(pointers,PyArray_INT,1,1);

   if (ip==NULL){
    PyErr_SetString(PyExc_ValueError, "Arg 2 did not make it to Int32 array");
    return NULL;
   }

   x=(PyArrayObject *)
      PyArray_ContiguousFromObject(vec,PyArray_FLOAT,1,1);

   if (matrix==NULL){
    PyErr_SetString(PyExc_ValueError, "Arg 3 did not make it to Float32 array");
    return NULL;
   }

   /* Check arguments, dim(IP)>=N and dim[A]>=IP[n] */

   if (ip->dimensions[0] < dim[0]) {
    PyErr_SetString(PyExc_ValueError, "Must have Arg 4 <= Arg2.shape[0]");
    return NULL;
   }
   if (x->dimensions[0] < dim[0]) {
    PyErr_SetString(PyExc_ValueError, "Must have Arg 4 <= Arg3.shape[0]");
    return NULL;
   }
   j= (* (int *) (ip->data + (dim[0]-1) * ip->strides[0]));

   if (matrix->dimensions[0] < j) {
    PyErr_SetString(PyExc_ValueError, "Arg1 not big enough for arg2[arg4]");
    return NULL;
   }
   result=(PyArrayObject *)PyArray_FromDims(1,dim,PyArray_FLOAT);
   if (result==NULL){
    PyErr_SetString(PyExc_ValueError, "Wow, could not make a results array");
   }
   solve_((REAL) matrix->data,
          (INTEGER) ip->data,
          (REAL) x->data,
          (REAL) result->data,
          dim);

   PyArray_XDECREF(matrix);
   Py_XDECREF(array);
   Py_XDECREF(vec);
   PyArray_XDECREF(ip);
   PyArray_XDECREF(result);
   PyArray_XDECREF(x);
   return PyArray_Return(result);
}




static PyMethodDef pylibcholmethods[] = {
   {"fcholdc", (PyCFunction) fcholdc, METH_VARARGS ,
   fcholdc_docstring},   
   {"fmatvec", (PyCFunction) fmatvec, METH_VARARGS ,
   fmatvec_docstring},
   {"fsolve",  (PyCFunction) fsolve , METH_VARARGS ,
   fsolve_docstring},
   {"getrefcount",  (PyCFunction) getrefcount , METH_VARARGS ,
   getrefcount_docstring},


   {NULL, NULL, 0, NULL} /* setinel */
};

void initpylibchol(void)
{
   PyObject *m, *d;
   m=Py_InitModule("pylibchol",pylibcholmethods);
   import_array();
   d=PyModule_GetDict(m);
}
@}




\section{Python wrapper}


A python wrapper for that stuff.

The following constructor just reads in a file. If you were to
inherit from this class you would just replace the constructor
with something else and get the rest of the functionality.

@D ciiconstructor
@{
   def __init__(self,filename):
      """
      Makes a python object representing the information coming from
      a DILS refinement in the MPRODD computer program
      Has a matrix which is the least squares matrix to be inverted
      during a Pawley refinement and also a set of non-intensity
      variables
      """
      f=open(filename,"r")
      items=f.readline().split()
      self.ni=int(items[0])
      self.nv=int(items[3])
      self.nt=self.ni+self.nv
      junk=f.readline()
      self.Ihkl=Numeric.zeros(self.ni,Numeric.Float32)
      self.HKL=Numeric.zeros((self.ni,3),Numeric.Float32)
      self.IP=Numeric.zeros(self.nt,Numeric.Int)
      for i in range(self.ni):
        try:
         items=f.readline().split()
         self.IP[i]=int(items[1])
         self.HKL[i]=[int(items[3]),int(items[4]),int(items[5])]
         self.Ihkl[i]=float(items[-1])/100
        except ValueError:
           print items
           print "Sorry - problems interpreting your file"
           sys.exit()
      self.rhs=Numeric.zeros(self.nt,Numeric.Float32)
      self.rhs[0:self.ni]=self.Ihkl
      for i in range(self.nv):
         items=f.readline().split()
         self.IP[self.ni+i]=int(items[1])
      print f.readline(),   
      self.matrix=Numeric.zeros(self.IP[-1],Numeric.Float32)
      self.matrix[0]=float(f.readline())
      for i in range(1,self.IP[-1]):
         self.matrix[i]=float(f.readline())
      f.close()
      print "Number of entries in matrix",self.IP[-1]
      print "Number of intensities=",self.ni,"Total vars",self.nt
      print "Percentage full=",self.IP[-1]*100.0/self.ni/self.ni
      self.ciimat=None
      self.L=None
      self.vecy=Numeric.zeros(self.nt,Numeric.Float32)
@}

This is the more useful part - it allows various things to be done
with the matrix.

@O cii.py
@{

@< pycopyright @>

import Numeric,math,sys
import pylibchol

class ciimatrix:
@< ciiconstructor @>

   def makesquare(self):
      """
      Produces a dense square matrix from the sparse one
      which is produced by the refinement program
      """
      mat=Numeric.zeros((self.nt,self.nt),Numeric.Float32,savespace=1)
      mat[0,0]=self.matrix[0]
      for i in range(1,self.nt):
        j=0
        while (self.IP[i]-j) > self.IP[i-1]:
           try:
              mat[i-j,i]=self.matrix[self.IP[i]-j-1]
              mat[i,i-j]=self.matrix[self.IP[i]-j-1]
              j+=1
           except IndexError:
              print i,j,self.IP[i],self.IP[i-1],\
                ciimat.shape,matrix.shape,self.IP[i]-j,i-j
              print "Trouble making a square matrix from the sparse one"
              sys.exit()
      return mat

   def diagonal(self):
      """
      Returns an array of the diagonal elements of the matrix
      """
      s=Numeric.zeros(self.nt,Numeric.Float32)
      s[0]=self.matrix[0]
      for i in range(self.nt): s[i]=self.matrix[self.IP[i]-1]
      return s


   def matvec(self,vec):
      """
      Python implementation of the sparse matrix vector computation
      FIXME - I think this has a bug?
      """
      print "Dont do this - it is really slow"
      res=Numeric.zeros(vec.shape[0],Numeric.Float32)
      res[0]=self.matrix[0]*vec[0]
      print res.shape
      for i in range(1,self.nt):
         # diagonal
         res[i]+=self.matrix[self.IP[i]-1]*vec[i]
         for j in range(self.IP[i-1],self.IP[i]-1):
            res[i]               +=vec[i-self.IP[i]+j+1]*self.matrix[j] # rows
            res[i-self.IP[i]+j+1]+=vec[i]        *self.matrix[j]        # columns
      return res    

   def fastmatvec(self,vec):
      """
      Fortran implementation of self.matvec
      """
      self.vecy=vec.astype(Numeric.Float32)
      return pylibchol.fmatvec(self.matrix,self.IP,self.vecy,self.nt)

   def solve(self,vec):
      """
      Uses the Cholesky decomposition, hopefully formed already to solve
      A x = b via 
      L (L^T x) = b    => L y = b solved for y
      then L^T x = y   => solved for x
      Uses self.L and returns x for a given b
      """
      if self.L==None: self.formchol()
      return pylibchol.fsolve(self.L,self.IP,vec,self.nt)
      
   def formchol(self,damp=None):
      """
      Forms Cholesky decomposition in self.L
      """
      if damp==None:
         m=self.matrix
      else:
         m=self.matrix.copy()
         for i in range(self.nt): m[self.IP[i]-1]*=(1.0+damp)
      self.L=pylibchol.fcholdc(m,self.IP,self.nt)

   def sparsechi2(self,vec):
      """
      Compute a correlated integrated intensity chi**2 for test vector vec
      chi**2 = (vec_obs - vec_calc)(M)(vec_obs - vec_calc)
      Uses the sparse matrix
      """
      omc=self.rhs-vec
      sum=Numeric.dot(omc,self.fastmatvec(omc))
#      print sum, sum/self.nt, math.sqrt(sum/self.nt)
      return sum

   def densechi2(self,vec,mat=None):
      """
      Compute a correlated integrated intensity chi**2 for test vector vec
      chi**2 = (vec_obs - vec_calc)(M)(vec_obs - vec_calc)
      Uses the full dense matrix
      """
      print "Dont do this - it is really slow!"
      omc=self.rhs-vec
      if mat==None: mat=self.makesquare()
      sum=Numeric.dot(omc,Numeric.dot(mat,omc))
#      print sum, sum/self.nt, math.sqrt(sum/self.nt)
      return sum
    
def normalise(matrix):
      """
      Does not really belong in this object.
      Makes diagonal elements equal to one and returns the scales
      (kills rows which cannot be normalised)
      """
      s=Numeric.zeros(matrix.shape[0],Numeric.Float32,savespace=1)
      for i in range(matrix.shape[0]):
         if (1.0+matrix[i,i]) > 1.0:
            s[i]=math.sqrt(matrix[i,i])
            matrix[:,i]=matrix[:,i]/s[i]
            matrix[i,:]=matrix[i,:]/s[i]
         else:
            s[i]=1.
            matrix[:,i]=0.
            matrix[i,:]=0.
            matrix[i,i]=1.
#         print "In Normalise",i,s[i]
      return matrix,s



if __name__=="__main__":
   import time
   t0=time.time()
   testobject=ciimatrix(sys.argv[1])
   t1=time.time()
   print "Time to read=",t1-t0
   print "\nSolving stuff\n"
   testobject.formchol(damp=0.1)
   t2=time.time()
   print "Time to form Cholesky decomposition=",t2-t1
   new=testobject.solve(test.fastmatvec(test.rhs)) # should give RHS back
   print test.sparsechi2(new),new[:10],"...\n\t\t",new[-3:]
   for i in xrange(10):
      new=new+testobject.solve(test.fastmatvec(test.rhs-new))
      print testobject.sparsechi2(new),new[:10],"...\n\t\t",new[-3:]
   omc=testobject.rhs-new

   aomc=test.fastmatvec(omc)
   for i in xrange(0,test.nt-5,5):
      print "i",i,Numeric.dot(omc[i:i+5],aomc[i:i+5])
   print omc[:10]
   print aomc[:10]
   sys.exit()
@}

\subsection{An optimisation object for correlated integrated intensities}

For refinement and plotting we only need to support the interface defined
by the powderdata object. Thus far that is only to contain x and y arrays 
used for plotting, a setrange function to decide what data to use, and 
a weightmatvec function for optmisation.

For computing models based on a set of correlated integrated intensities
we might need some extra functionality, which will be added here as needed.

@o ciidata.py
@{
""" Correlated integrated intensity object for use in refinement """
@< pycopyright @>

from Numeric import *
import cii,pylibchol

class ciidata:
   def __init__(self,ciiobject ,d):
      self.ciiobject=ciiobject   # A cii object as above
      self.d=d                   # Holds metadata in dictionary
      self.x=self.xfull=array(range(self.ciiobject.ni)) 
                     # default to number of integrated intensities
      self.y=self.yfull=self.ciiobject.Ihkl
      self.npoints=self.ciiobject.ni
      try: 
        self.rcm = d['rcm']
      except KeyError:
        print "No unit cell supplied for your data"
      self.lh=None
   def setrange(self,lh=None): 
      self.lh=lh
      if self.lh==None:
         self.x=self.xfull
         self.y=self.yfull
         self.npoints=self.x.shape[0]
      else:         
         l=0 # always - not sure how to throw out low angle peaks yet
         h=searchsorted(self.xfull,self.lh[1])
         self.x=self.xfull[l:h]
         self.y=self.yfull[l:h]
         self.d["active_range"]=[l,h]
         self.lh=[l,h]  
         self.npoints=self.x.shape[0]
   def weightmatvec(self,vector):
      v=vector.astype(Float32)
      return pylibchol.fmatvec(self.ciiobject.matrix,self.ciiobject.IP,v,self.npoints)
   def sinthetaoverlambda2(self,i):
      hkl=self.ciiobject.HKL[i]
      return dot(hkl,matrixmultiply(self.d['rcm'],hkl))
   def setfsq(self,fsq):
      self.fsq=fsq # holds list of computed structure factors if available

@}
  

Now for something useful - we want to refine solvent scattering (later on)
and need to read a dils file from prodd, a ccl file to get the cell parameters
and a pdbfile to get the structure factors.

@o solventrefinedata.py
@{
""" Data for refining a bulk solvent model """
@< pycopyright @>
from Numeric import *
import ciidata,cii
from math import cos,pi
from LinearAlgebra import inverse   
import cctbxfcalc  # For a set of fcalc numbers

class solventrefinedata(ciidata.ciidata):
   def __init__(self,dilsfile, cclfile, pdbfile):
      ciiobj = cii.ciimatrix(dilsfile)
      fcalc = cctbxfcalc.getfcalcfrompdb(pdbfile,ciiobj.HKL.astype(Int).tolist())
      fsq=abs(fcalc).astype(Float32)
      fsq=fsq*fsq*1e-6
      d={}
      c=open(cclfile,"r").readlines()
      for line in c:
         if line[0]=="C":
            [a,b,c,alpha,beta,gamma]=map(float,line[1:].split())
      d['rm'] = array( [ [ a*a, a*b*cos(gamma*pi/180) , a*c*cos(beta *pi/180) ] ,
                         [ a*b*cos(gamma*pi/180) , b*b, b*c*cos(alpha*pi/180) ] ,
                         [ a*c*cos(beta *pi/180) , b*c*cos(alpha*pi/180) , c*c] ] )
      d['rcm']= inverse(d['rm'])
      ciidata.ciidata.__init__(self,ciiobj ,d)
      self.setfsq(fsq)
      
@}


\chapter{Building models}

Allow the user to pick the crystal structure, peak shape, correction
factors and so on and so forth and put them together in any 
way they choose.
This is a generalised way to chain together pieces of mathematics
to make a full model.
Here for deciding which crystallographic phases are present in 
which histograms with which kinds of peakshapes and absorbtion
corrections. etc.

The models will be put together using sums and multiplications.
An example would be that the data are the background plus a 
scale factor times a correction factor times a crystal phase.
Derivatives are then added or adjusted using the chain rule.
We can use a tree structure to represent this model.
\[ y_c = b + scI \]
\[ \frac{\partial y_c}{\partial x} = \frac{\partial b}{\partial x} 
        +cI\frac{\partial s}{\partial x}
        +sI\frac{\partial c}{\partial x} 
        +sc\frac{\partial I}{\partial x} \]
If we know which parameters have derivatives with respect to which
variables, then certain terms in the summation above will disappear.
So each contributor to a model needs to uniquely identify a list
of variables which it is using.

A contributor then offers a value and the derivatives of the value
with respect to a set of parameters which it tells you about.
A model is the sum of a set of contributors, and contributors may 
also be multiplied together.
A crystal phase might be seen as a sum of contributors (the peaks)
where the derivatives of intensities might be all with respect to 
different parameters (pattern decompositions) or all with respect
to the same set of parameters (crystal structure refinement).

FIXME - read and understand the parts of Programming python on trees
before going on. We need to find a way to do the computations efficiently
without sacrificing generality.

\section{Models viewed from the optimisation perspective}

The model will deal with all the crap about parameters and error propagation
from refined things into varied things. 
A model object is only going to supply the computed value for a datapoint
and derivatives of that value for a set of variables.
In order to cope with sparse models in a simple way it will supply a 
list of derivatives and a list of pointers indicating which 
derivative is which.

For example, fitting a load of peaks on a background function, most points
only have derivatives for one or two peaks, but all points have derivatives
for background. The model object has a memmber indicating the number
of variables in total.
When it is called for a particular ordinal (x) value it returns with 
$y_{calc}$ and derivatives $dy_{calc}/dx_i$ for the set of variables
$x_i$ having non-zero derivatives at this point. 
An integer array is needed giving the list of $i$ values so we know
which derivative is which.

For certain peak shape functions it may be more convenient to compute the
entire peak and store it in memory so that the model is not called for 
each point once at a time, but for an array of ordinal points. 
In that case it needs to return the full list of derivatives for each point. 
In practice this could be wrapped by making the full list on the first 
call and then returning what is needed, or having a rolling buffer in the
model.

For observations which are correlated several datapoints are needed at the 
same time, but weights are to be treated by the optimisation object, so
which derivatives are required are determined by the optimiser deciding
which x-points it calls the model for.


A simple class to implement a model is then:

@o model.py
@{
@< pycopyright @>
from Numeric import *
class model:
   """
   Intended as a base class for inheriting - you need to set up a list
   called self.variables and it will set up arrays to hold them and 
   handle most of the talking to optimisers
   """
   def __init__(self,**kwds):
      self.vd={}         # Dictionary to find variables in variable value array (vv)
      self.set={}        # Dictionary indicating whether variables are initialised
      i=0
      for item in self.variables: 
         self.vd[item]=i        # Number the variables into arrays
         self.set[item]=False   # Indicate values have not been initialised
         i+=1
      # variable values and error bar array
      self.vv=zeros(len(self.variables),Float32)
      self.ve=zeros(len(self.variables),Float32)
      # Fill out any variable values which are present
      for key in kwds.keys():
         try: 
            self.vv[self.vd[key]] = kwds[key]
            self.set[key]=True
         except KeyError:
            print "Unrecognised Keyword argument:",key,"=",kwds[key]
      self.ycalc=None

   def compute(self,data):
      pass

   def gradient(self,variable):
      pass

   def apply_shift(self,variable,shift):
      self.vv[self.vd[variable]] = self.vv[self.vd[variable]] + shift

   def set_value(self,variable,value):
      self.vv[self.vd[variable]] = value
      self.set[variable]=True

   def get_value(self,variable):
      return self.vv[self.vd[variable]]

   def set_errorbar(self,variable,value):
      self.ve[self.vd[variable]] = value

   def get_errorbar(self,variable):
      return self.ve[self.vd[variable]] 


   def get_variables(self):
       return self.variables
@}

Then to refine a single constant value we only need to fill out the list
of variables being used and override the compute and gradient methods.

@o constantmodel.py
@{
@< pycopyright @>
import modelclass
from Numeric import *
class constantmodel(model):
   def __init__(self,**kwds):
      """ A single constant value, call constantmodel(constant=y) """
      self.variables["constant"]
      model.__init__(**kwds)

   def compute(self,data):
      self.ycalc=ones(data.y.shape,Float32)*self.vv[self.vd['constant']]

   def gradient(self,variable):
      return ones(self.ycalc.shape,Float32)
@}

etc...

\section{Slow but flexible models}

Given the model class of the previous section, we can make
another model class which calls on the one above. 
For example, if we want to define a function to be fitted as
a sum of two things we can make a "summing" model class. This
is initialised with the two things it is summing and 
passes the derivatives on appropriately.

eg:
@o summodel.py
@{
@< pycopyright @>
class summodel(model):
   def __init__(m1,m2):
      """ Sums the two models m1 and m2 together """
      self.nv = m1.nv + m2.nv # assume independent???
      self.m1=m1
      self.m2=m2
   def compute(self,x):
      m1yc = self.m1.compute(x) 
      m2yc = self.m2.compute(x)
      yc = m1yc[0]+m2yc[0]
      iw = self.join(m1yc[1],m2yc[1]) 
      di = self.join(m1y2[2],m2yc[2])
   def applyshift(self,shift):
      ???
@}

Need to make unique lists of variables somehow.


\chapter{Optimisation}

A chapter going into the specialised details of how to do the 
optimisation with the interfaces defined in the data structures
chapter.
Eventually develop c or fortran routines which can do the optimisation
efficiently.
Broadly two routines are used. 
One adds in ``observations'' to the least squares refinement, the 
other finds a solution vector to the least squares problems and 
determines some sort of error and covariance estimates.

For refinements which go through a number of cycles we should
use the covariance information to find linear combinations of parameters
which are not so correlated with each other. (Automatic refinement using
eigenvectors of variables - should go back to find out if these are
not approximately the basic variables to take sum and difference
earlier on and avoid rounding errors).



\section{Dense gradients, 1D or 2D weight matrix}

This started as the example in the design section, with plans for
a fortran routine which worried about the matrix being symmetric
and triangular.

After some fooling around fitting peaks to 1D datasets this
ends up being an optimiser when all gradients are present
at the same time for all data points (ie the scratch file
in prodd)

This merges in the planned 2D weight matrix optmiser by calling on 
the data object to possess a member function which offers
a weighted-matrix-vector product. For the 1D stuff that it is 
trivial, and in some ways more convenient to move the treatment
of statistical weights into the data object.


@o densemodelfit.py
@{
""" Least squares optimisation for dense models and 1D weights """
@< pycopyright @>
import time
from Numeric import *
from LinearAlgebra import generalized_inverse
class densemodelfit:
   def __init__(self,model=None,data=None,marq=None):
      """
      Optimises a model against a set of data
      The data should provide array members data.y and data.e which
        are used for the observed data values and error bars
      The model should generate an array member model.ycalc of the
        same dimension as the data.y. It might call on other information
        in the data object to do that.
        The model must also supply a list of variables (strings) accessible
        via a "get_variables()" member function.
        It must also supply a function to apply shifts
        Mainly the model supplies a "compute" member function which creates
        a member ycalc array (the computed model) and a 2D array of gradients
        holding d{model.ycalc}/d{variable}. Positions of the variables in
        the gradients array is given by the member model.grad_loc dictionary.
      """
      self.data=data
      self.model=model
      self.variables=model.get_variables() # make a copy
      self.vd = {} # dictionary for labelling variables
      i=0
      # Positions of variables in least squares matrix
      for item in self.variables: self.vd[item]=i; i+=1
      self.nv=len(self.variables)
      self.icyc=0
      if marq==None: 
         self.marq=1.0
      else:
         self.marq=marq

   def refine(self, ncycles=1):
      start=time.time()
      end=self.icyc+ncycles
      while self.icyc < end: 
         self.model.compute(self.data)   # model only knows about data when computing
         obsminuscalc=self.data.y-self.model.ycalc
# 1D         ge=zeros((self.data.npoints,self.nv),Float32)               # gradient array
# 1D        for i in range(self.nv):
# 1D              ge[:,i] = (self.model.gradient(self.variables[i]) / self.data.e).astype(Float32)
# 1D        self.lsqmat=matrixmultiply(transpose(ge),ge)
# 1D        self.rhs=matrixmultiply(obsminuscalc/self.data.e,ge)
         g=zeros((self.data.npoints,self.nv),Float32)               # gradient array
         for i in range(self.nv):                        # avoid this somehow?
            g[:,i] = self.model.gradient(self.variables[i])
         wg=zeros((self.data.npoints,self.nv),Float32)     # weighted gradient array
         for i in range(self.nv):
            wg[:,i] = self.data.weightmatvec(g[:,i])
         self.rhs=matrixmultiply(obsminuscalc,wg)   
         self.lsqmat=matrixmultiply(transpose(g).astype(Float),wg.astype(Float))
         # Error bars
         self.emat=self.generalized_inverse(self.lsqmat)
         # damping:
         if self.marq != 1.:
            for i in range(self.nv):
               self.lsqmat[i,i] = self.marq*self.lsqmat[i,i]
         # invert to find shift (no need really...)
         self.inverse=self.generalized_inverse(self.lsqmat)
         shifts = matrixmultiply(self.inverse,self.rhs).astype(Float32)
# 1D          self.chi2=sum(obsminuscalc*obsminuscalc/self.data.e/self.data.e)
         self.chi2=sum(obsminuscalc* self.data.weightmatvec(obsminuscalc) )
         if self.data.npoints-self.nv > 0:
         self.reducedchi2=self.chi2/(self.data.npoints-self.nv)
            self.emat=self.emat*self.reducedchi2
         else:
            self.reducedchi2=0. # No chi^2 correction here - it makes errors wrongly zero
         print "Cycle ",self.icyc," chi^2 = ",self.chi2, "Reduced chi^2=",self.reducedchi2
         for item in self.variables:
            self.model.apply_shift(item,shifts[self.vd[item]])
         self.icyc+=1
      #
      # End cycles loop
      print "Variable  Value           Esd             Shift           Shift/esd"
      for item in self.variables:
         i=self.vd[item]
         e=sqrt(self.emat[i,i])
         self.model.set_errorbar(item,e)
         print "%8s  %-10.7e"%( item,self.model.get_value(item)),
         if e>0.:
            print " %-10.7e  %-10.7e  %-10.7e"%(e, shifts[i] , shifts[i]/e )
         else:
            print " %-10.7e  %-10.7e  NaN"%(e, shifts[i] )
#      print "Lsq Matrix", self.lsqmat
      print "Time for that cycle was:",time.time()-start
         
   def generalized_inverse(self,matrix):
       # marker to try out or not scaling the matrix
       # I hope the LinearAlgebra package does not need this??
       # The idea is to make the diagonal of the matrix have equal elements
       # which is equivalent to making all variables be in units of ~1 esd
       try:
       s=sqrt(diagonal(matrix))
       except ValueError:
          print "Least squares matrices should always (by definition) have positive",\
                "diagonal elements"
          print matrix
          print diagonal(matrix)
          raise
       mycopy=matrix.copy().astype(Float)
       for i in range(matrix.shape[0]):
          if s[i]>0.:
             mycopy[i,:]=mycopy[i,:]/s[i]
             mycopy[:,i]=mycopy[:,i]/s[i]
          else:
             mycopy[:,i]=0. ; mycopy[i,:]=0.; mycopy[i,i]=1. ; s[i]=1.
       inv = generalized_inverse(mycopy)   
       for i in range(matrix.shape[0]):
          inv[i,:]=inv[i,:]/s[i]
          inv[:,i]=inv[:,i]/s[i]
       return inv

@}

\section{Sparse gradients, 1D weight matrix}


This optimiser can extend the one in the design section to work 
with sparse gradients. There is a possible optimisation which is 
still needed to decide which matrix elements are needed in the
refinement matrix. 
If the gradients are sparse we need to find out which pairs 
of gradients need to have a least squares matrix element. 
This means building the sparse matrix structure which 
is used in mprodd - or perhaps a more general sparse matrix
structure with arbitrary off diagonal elements, with a 
fairly quick matrix vector product to use iterative
methods for solving the system (eg CGLS).

The optimiser will expect to be able to call on a model
object for the following things:
\begin{itemize}
\item Number of variables = model.nv
\item Computed function = model.compute(x) ; This will return a number $y_calc$,
and arrays iw and di
\item Apply shifts = model.applyshifts(shift)
\end{itemize}


\section{Dense gradients, 2D weight matrix}

This would be used for the optimisation of a crystal structure against
correlated integrated intensity data.

All that needs to be modified from the 1 dimensional weight matrix follows
from the mathematics described in the opening chapter. In the original
implementation of refinement against a Lanthanum Hexaboride structure
the key lines were as follows:

\begin{verbatim}
  for i in range(nv):
     rhs[i]=Numeric.dot(diff,myciidata.fastmatvec(g[:,i]))
     for j in range(nv):
         lsqA[i,j]+=Numeric.dot(myciidata.fastmatvec(g[:,i]),g[:,j])
         lsqA[j,i]+=Numeric.dot(g[:,j], myciidata.fastmatvec(g[:,i]))
  chi=Numeric.dot(diff,myciidata.fastmatvec(diff))
  rchi=chi/(nt-nv)
\end{verbatim}

When the weight matrix is diagonal this reduces to:

\begin{verbatim}
         for i in range(self.nv):
               g[:,i] = (self.model.gradient(self.variables[i]) / self.data.e).astype(Float32)
         self.lsqmat=matrixmultiply(transpose(g),g)
         self.rhs=matrixmultiply(obsminuscalc/self.data.e,g)
\end{verbatim}

In the implementation of the dense gradients and 1D data weight matrix we 
already have all the information needed, as it currently requires
derivatives of all data points to be present with respect to all parameters.
So we only need to modify the densemodels object sufficiently to take
the error bars into account properly. Indeed, we might note that a diagonal
weight matrix is merely a special case of a non-diagonal weight matrix
and combine both into the same method. There ought to be no trade-off 
involved in doing so?

So we need to transform the four lines in the fitdensemodels1d 
from being \code{g = model.gradient/data.e ; lsqmat=g x $g{^T}$ ; rhs=g.diff/e}
to something else. If we call a data.weightmatvec(vector) which returns
the vector multiplied by the weights, then we are away. 

If you are reading this document it will make little sense, but I leave this
section in with the answer to the apparent problem now implemented in the
densemodelfit object.


\begin{verbatim}
         g=zeros((self.data.npoints,self.nv),Float32)               # gradient array
         for i in range(self.nv):                        # avoid this somehow?
            g[:,i] = self.model.gradient(self.variables[i])
         wg=zeros((self.data.npoints,self.nv),Float32)     # weighted gradient array
         for i in range(self.nv):
            wg[:,i] = self.data.weightmatvec(self.model.gradient(self.variables[i]))
         self.rhs=matrixmultiply(obsminuscalc,wg)   
         self.lsqmat=matrixmultiply(transpose(g),wg)
\end{verbatim}

And with that the author was able to shut down the laptop before midnight with
a satisfied smile. The correlations have been tamed.


We can check that this is correct (or not) by setting up a test model with
a test dataset.

@o testdensedemodelfit.py
@{
""" Testing for the non-diagonal weight matrix least squares """
@< pycopyright @>

import model, densemodelfit
from Numeric import *
class testmodelclass(model.model):
   def __init__(self,**kwds):
      self.variables = map(str, range(kwds['npoints']) )         
      model.model.__init__(self)
   def compute(self,data):
      self.ycalc=self.vv
      self.g=zeros(self.vv.shape,Float32)
   def gradient(self,variable):
      g=zeros(self.vv.shape,Float32)
      g[self.vd[variable]]=1.
      return g

class testdataclass:
   def __init__(self):  
      import RandomArray
      self.x=arange(3)
      self.y=RandomArray.random(3) # obs
      self.npoints=3
      self.matrix = array( [[ 1.0, 1.0, 0.1] , [ 1.0, 1.0, 0.1] , [0.1, 0.1, 2.0 ] ], Float32)
      self.d={}
   def weightmatvec(self,vec):
      return matrixmultiply(self.matrix,vec)

if __name__=="__main__":
   m=testmodelclass(npoints=3)
   d=testdataclass()
   o=densemodelfit.densemodelfit(data=d,model=m)
   o.refine(ncycles=3)
   print "Data :",d.y
   print "Model:",m.ycalc
   diff=m.ycalc-d.y
   print "Uncorrelated residual",sum(diff*diff)
@}

\section{Sparse gradients, 2D weight matrix}

This would be used for optimisation against a set of correlated integrated
intensities, but the gradients can become sparse for example in a multiphase
refinement, where one phase is protein and the other is not.



\chapter{Peakshapes}

Wide variety of peakshapes exist.
Look out the python wrapper to Larry Finger's low angle asymmetry
correcting routine for neutron, synchrotron and lab x-rays with
a constant wavelength.
Make a wrapper for Bill David's CCSL ToF peakshape, or implement
the peakshapes in GSAS.

Aim to be slightly original here and make the peak broadening
effects as a function of scattering vector (Q) to wrap onto
whatever the peakshape function. 
Should make things the same whether using ToF, CW or whatever.
Only problem is to describe instrumental broadenings in Q, which 
may be very inconvenient? Look into it.

\section{Constant wavelength}

\subsection{Gaussian peakshape}

A simple peakshape function for comparison and debugging purposes
mainly. Implement a Gausssian here FIXME

\subsection{Larry Finger's low angle asymmetry correcting subroutine}

As a quick fix we will put Larry Finger's asymmetry correcting
routine so that we have something available to get started.
A wrapper was written before, which is included here.

\subsubsection{The python-c wrapper}

The c-wrapper simply calls the function passing args forward
and returning result values. Since the arguments to the function are
either in or out going, it only takes the input ones and only returns
the output ones. The prototype for the function assumes g77/gcc link
conventions. Since the author has access to those compilers on both
windows and linux that seems like a nice start. One day it 
might be worth finding out about the new generation of unix
based macintosh computers.

A second version could be made which takes a numeric array of twotheta
values to compute the function on, and which could returns arrays
of calculated values and derivatives... FIXME


@O cprofval.c
@{
@< ccopyright @>

#include <Python.h>                  /* To talk to python */
#include "Numeric/arrayobject.h"

static char docstring[] =
"Larry Fingers low angle correcting peakshape\n"\
"profval:\n"\
"Call with float args: Eta,Gamma,S_L,D_L,TwoTH,TwoTH0,use_asym[int: 0/1]\n"\
"Returns function_value, dPRdT, dPRdG, dPRdE , dPRdS , dPRdD, area\n"\
"\n"\
"profval_array: As above but call with arrya of TwoTH values and \n"\
"get back arrays of function and derivative values\n"\
"\n"\
"! fortran comments follow\n"\
"       real*4 function Profval( Eta , Gamma , S_L , D_L , TwoTH , \n"\
"     1   TwoTH0 , dPRdT, dPRdG, dPRdE , dPRdS , dPRdD , Use_Asym )\n"\
"c\n"\
"c Returns value of Profile\n"\
"c   Eta is the mixing coefficient between Gaussian and Lorentzian\n"\
"c   Gamma is the FWHM\n"\
"c   S_L is source width/detector distance\n"\
"c   D_L is detector width/detector distance\n"\
"c   TwoTH is point at which to evaluate the profile\n"\
"c   TwoTH0 is two theta value for peak\n"\
"c   dPRdT is derivative of profile wrt TwoTH0\n"\
"c   dPRdG is derivative of profile wrt Gamma\n"\
"c   dPRdE is derivative of profile wrt Eta\n"\
"c   dPRdS is derivative of profile wrt S_L\n"\
"c   dPRdD is derivative of profile wrt D_L\n"\
"c   Use_Asym is true if asymmetry to be used\n"\
"c\n"\
"c\n"\
"c Asymmetry due to axial divergence using the method of Finger, Cox and\n"\
"c    Jephcoat, J. Appl. Cryst. 27, 892, 1992.\n"\
"\n"\
"      real*4 Eta , Gamma , S_L , D_L , TwoTH \n"\
"      real*4 TwoTH0 , dPRdT, dPRdG, dPRdE , dPRdS , dPRdD\n"\
"      logical Use_Asym\n";

/* c prototype */
float profval_(float *eta,   float *gamma,  float *s_l,   float *d_l, 
               float *twoth, float *twoth0, float *dprdt, float *dprdg, 
	            float *dprde, float *dprds,  float *dprdd, int *use_asym);

/* to be called from python */
static PyObject * profval (PyObject *self, PyObject *args)
{
   float eta, gamma, s_l, d_l, twoth, twoth0, area; /* input */
   int use_asym; /* are int and logical always compatible? */
   
   float result, dprdt, dprdg, dprde, dprds, dprdd;   /* output */
   
   if(!PyArg_ParseTuple(args, "fffffffi", &eta, &gamma, &s_l, &d_l, 
                        &twoth, &twoth0, &use_asym, &area) )   return NULL;
			
   result=profval_(&eta, &gamma, &s_l, &d_l, &twoth, &twoth0, 
                   &dprdt, &dprdg, &dprde, &dprds, &dprdd, &use_asym);

   return Py_BuildValue("ffffff",result*area, 
                                 dprdt*area, dprdg*area, dprde*area, dprds*area, dprdd*area);
}   

static PyObject * profval_array (PyObject *self, PyObject *args)
{
/* Exactly the same as above, but called on array of TwoTH values
   and returns arrays of derivatives */
   float eta, gamma, s_l, d_l, twoth, twoth0; /* input */
   int use_asym,i,npts,dim[1]; /* are int and logical always compatible? */
   PyArrayObject *ttharray, *result_a, *dprdt_a, *dprdg_a, *dprde_a, *dprds_a, *dprdd_a;
   float result, dprdt, dprdg, dprde, dprds, dprdd, area;   /* output */
   
   if(!PyArg_ParseTuple(args, "ffffO!ffi", &eta, &gamma, &s_l, &d_l, 
                        &PyArray_Type, &ttharray,
                        &twoth0, &area,&use_asym) )   return NULL;
   npts=ttharray->dimensions[0] ;
   dim[0]=npts;
   /* Check type of ttharray is OK */
   if(ttharray->descr->type_num != PyArray_FLOAT ){
      PyErr_SetString(PyExc_ValueError,"Two theta array must be 1D Numeric.Float32 type") ;
                        return NULL;
   }
   if(ttharray->nd != 1 ){
      PyErr_SetString(PyExc_ValueError,"Two theta array must be 1D Numeric.Float32 type") ;
                        return NULL;
   }


   /* Make results arrays */
   result_a = (PyArrayObject *)PyArray_FromDims(1,dim,PyArray_FLOAT);
   if(result_a == NULL){PyErr_SetString(PyExc_ValueError,"Couldnt allocate for results") ;
                        return NULL; }
   dprdt_a = (PyArrayObject *)PyArray_FromDims(1,dim,PyArray_FLOAT);
   if(dprdt_a == NULL){PyErr_SetString(PyExc_ValueError,"Couldnt allocate for results") ;
                        return NULL; }
   dprdg_a = (PyArrayObject *)PyArray_FromDims(1,dim,PyArray_FLOAT);
   if(dprdg_a == NULL){PyErr_SetString(PyExc_ValueError,"Couldnt allocate for results") ;
                        return NULL; }
   dprde_a = (PyArrayObject *)PyArray_FromDims(1,dim,PyArray_FLOAT);
   if(dprde_a == NULL){PyErr_SetString(PyExc_ValueError,"Couldnt allocate for results") ;
                        return NULL; }
   dprds_a = (PyArrayObject *)PyArray_FromDims(1,dim,PyArray_FLOAT);
   if(dprds_a == NULL){PyErr_SetString(PyExc_ValueError,"Couldnt allocate for results") ;
                        return NULL; }
   dprdd_a = (PyArrayObject *)PyArray_FromDims(1,dim,PyArray_FLOAT);
   if(dprdd_a == NULL){PyErr_SetString(PyExc_ValueError,"Couldnt allocate for results") ;
                        return NULL; }

   /* Loop over the data - all this crap to put a loop into c! */
   for(i=0; i<npts ; i++){
	   twoth=(* (float *) (ttharray->data + i*ttharray->strides[0]) );
      result=profval_(&eta, &gamma, &s_l, &d_l, &twoth, &twoth0, 
                      &dprdt, &dprdg, &dprde, &dprds, &dprdd, &use_asym);
      (* (float *) (result_a->data + i*result_a->strides[0])) = result*area;
      (* (float *) (dprdt_a->data  + i* dprdt_a->strides[0])) = dprdt*area;
      (* (float *) (dprdg_a->data  + i* dprdg_a->strides[0])) = dprdg*area;
      (* (float *) (dprde_a->data  + i* dprde_a->strides[0])) = dprde*area;
      (* (float *) (dprds_a->data  + i* dprds_a->strides[0])) = dprds*area;
      (* (float *) (dprdd_a->data  + i* dprdd_a->strides[0])) = dprdd*area;
      }
   /* Decref everything we allocated here as we won't look at it again */
   PyArray_XDECREF(ttharray);
   PyArray_XDECREF(result_a);
   PyArray_XDECREF(dprdt_a);
   PyArray_XDECREF(dprdg_a);
   PyArray_XDECREF(dprde_a);
   PyArray_XDECREF(dprds_a);
   PyArray_XDECREF(dprdd_a);

   /* Use "N" in Py_BuildValue to return the results. "O" increfs */
   /* returning a tuple of arrays */
   return Py_BuildValue("NNNNNN",  PyArray_Return(result_a),PyArray_Return(dprdt_a),
          PyArray_Return(dprdg_a), PyArray_Return(dprde_a), PyArray_Return(dprds_a), 
          PyArray_Return(dprdd_a) );
}   


/* For making the python extension */
static PyMethodDef profvalMethods[] = {
   {"profval", (PyCFunction) profval, METH_VARARGS, docstring},
   {"profval_array", (PyCFunction) profval_array, METH_VARARGS, docstring},
   {NULL, NULL, 0, NULL} /* setinel */
};

/* More gubbins for making the python extension */
void initprofval(void)
{
   PyObject *m, *d;
   m=Py_InitModule("profval", profvalMethods);
   import_array()
   d=PyModule_GetDict(m);
}
@}

\subsubsection{The fortran code}

The routine has been altered to return derivatives of the 
sum and difference of the asymmetry parameters. 
The large table of precomputed parameters for doing numerical
integrations are placed in the appendix to make the
text less painful to read.

@O profval.f
@{        real*4 function Profval( Eta , Gamma , S_L , D_L , TwoTH , 
     1   TwoTH0 , dPRdT, dPRdG, dPRdE , dPRdS , dPRdD , Use_Asym )
c  [119.19 bugfixed? jpw 15-oct-2001, Feb 02 yet another thing]
c
c Thanks to LWF for passing on those fixes and assistance with
c using the routine, and most of all, for making it available!
c
c 
c
c Returns value of Profile
c   Eta is the mixing coefficient between Gaussian and Lorentzian
c   Gamma is the FWHM
c   S_L is source width/detector distance
c   D_L is detector width/detector distance
c   TwoTH is point at which to evaluate the profile
c   TwoTH0 is two theta value for peak
c   dPRdT is derivative of profile wrt TwoTH0
c   dPRdG is derivative of profile wrt Gamma
c   dPRdE is derivative of profile wrt Eta
c   dPRdS is derivative of profile wrt S_L
c   dPRdD is derivative of profile wrt D_L
c   Use_Asym is true if asymmetry to be used
c
c Asymmetry due to axial divergence using the method of Finger, Cox and
c    Jephcoat, J. Appl. Cryst. 27, 892, 1992.
      implicit none
      real*4 Eta , Gamma , S_L , D_L , TwoTH 
      real*4 TwoTH0 , dPRdT, dPRdG, dPRdE , dPRdS , dPRdD
      logical Use_Asym
      integer*4 NTERMS(14)/6,10,20,40,60,80,100,150,200,300,400,
     1   600,800,1000/
C 119.19 changed first term to 1 instead of 0 - JPW
      integer*4 Fstterm(14)/1,3,8,18,38,68,108,158,233,333,483,
     1   683,983,1383/
      real*4 RAD/57.2957795/
      integer*4 ArrayNum , K , NGT, ngt2 , it, i
      real*4 CsTH             	! cos(theta)
      real*4 TTH		! tan(theta)
      real*4 SnTwoTH		! sin(twoth)
      real*4 CsTwoTH 		! cos(twoth)
      real*4 ApB		! (S + H)/L 
      real*4 AmB		! (S - H)/L 
      real*4 ApB2 		! (ApB) **2
      real*4 Einfl              ! 2phi value for inflection point 
      real*4 Emin               ! 2phi value for minimum 
      real*4 dEmindA            ! derivative of Emin wrt A
      real*4 tmp , tmp1 , tmp2  ! intermediate values 
      real*4 WP(1883) , XP(1883)! Storage for Gauss-Legendre weights and intervals 
      real*4 Delta              ! Angle of integration for comvolution 
      real*4 dDELTAdA           ! derivative of DELTA wrt A (S/L)
      real*4 sinDELTA           ! sine of DELTA 
      real*4 cosDELTA           ! cosine of DELTA 
      real*4 tanDELTA           ! tangent of DELTA 
      real*4 RcosDELTA          ! 1/cos(DELTA) 
      real*4 F , dFdA, G , dGdA , dGdB , PsVoigt,  stepsize 
      real*4 sumWG , sumWRG , sumWdGdA , sumWRdGdA ,sumWdGdB , sumWRdGdB
      real*4 sumWGdRdG , sumWGdRdE , sumWGdRdA , sumWGdRdB , sumWGdRd2t 
@<profvaldata@>
      CsTH = cos(TwoTH0 * 0.5/RAD)
      if (abs(CsTH) .lt. 1.0e-15) CsTH = 1.0e-15
      TTH = sin(TwoTH0 * 0.5/RAD)/CsTH
      CsTwoTH = cos(TwoTH0/RAD)
      SnTwoTH = sin(TwoTH0/RAD)
      ApB = S_L + D_L
      AmB = S_L - D_L
      ApB2 = ApB**2
      if (((S_L .ne. 0.0) .or. (D_L .ne. 0.0)) .and. Use_Asym) then
        tmp = sqrt(1.0 + AmB**2)*CsTwoTH
        if (abs(tmp) .gt. 1.0) then
          Einfl = acos(CsTwoTH)*RAD
        else
          Einfl = acos(tmp)*RAD
        endif
        tmp2 = 1.0 + ApB2
        tmp = sqrt(tmp2 ) * CsTwoTH

c If S_L or D_L are zero, set Einfl = 2theta 
    
        if ((S_L .eq. 0.0) .or. (D_L .eq. 0.0)) Einfl = TwoTH0
        if (abs(tmp) .le. 1.0) then
          Emin = acos(tmp) * RAD
          tmp1 = tmp2 * (1.0 - tmp2 * CsTwoTH**2)
        else
          tmp1 = 0.0
          if (tmp .gt. 0.0) then
            Emin = 0.0
          else
            Emin = 180.0
          endif
        endif
        if ((tmp1 .gt. 0.0) .and. (abs(tmp) .le. 1.0)) then
          dEmindA = -ApB * CsTwoTH/sqrt(tmp1)
        else
          dEmindA = 0.0
        endif
        ArrayNum = 1
        K = 400.0 * (TwoTH0 - Emin)   ! Calculate number of terms needed 
C From LWF - number of terms must be such that the interval between 2phi(min) and 2theta
C is in steps no larger than 0.005
C        stepsize = (twoth0-emin)/float(K)
        stepsize = 1.0/400.0 ! This seems to be always the case apart from rounding K
        if(stepsize .gt. gamma/10.0) then
           stepsize = gamma/10.0
           K = (twoth0-emin)/stepsize
        endif
        do while ((ArrayNum .lt. 14) .and. (K .gt. NTERMS(ArrayNum)))
          ArrayNum = ArrayNum + 1
        enddo
        NGT = nterms(ArrayNum)              ! Save number of terms 
        ngt2 = ngt / 2
c Clear terms needed for summations 
        sumWG = 0.0
        sumWRG = 0.0
        sumWdGdA = 0.0
        sumWRdGdA = 0.0
        sumWdGdB = 0.0
        sumWRdGdB = 0.0
        sumWGdRd2t = 0.0
        sumWGdRdG = 0.0
        sumWGdRdE = 0.0
        sumWGdRdA = 0.0
        sumWGdRdB = 0.0
c Compute the convolution integral 
        it = fstterm(arraynum)-ngt2
        do K = ngt2 , NGT 
          delta = Emin + (TwoTH0 - Emin) * xp(K + it)
          dDeltadA = (1.0 - xp(k+it) ) * dEmindA
          sinDELTA = sin(Delta/RAD)
          cosDELTA = cos(Delta/RAD)
          if (abs(cosDELTA) .lt. 1.0e-15) cosDELTA = 1.0e-15
          RcosDELTA = 1.0 / cosDELTA
          tanDELTA = tan(Delta/RAD)
          tmp = cosDELTA**2 - CsTwoTH**2
          if (tmp .gt. 0.0) then
            tmp1 = sqrt(tmp)
            F = abs(CsTwoTH) / tmp1
            dFdA = cosDELTA * CsTwoTH * sinDELTA * dDELTAdA 
     1         / (tmp1 * tmp1 * tmp1)
          else
            F = 0.0
            dFdA = 0.0
          endif
c  calculate G(Delta,2theta) , FCJ eq. 7a and 7b 
          if ( abs(Delta - Emin) .gt. abs(Einfl - Emin)) then
            if (S_L .gt. D_L) then
!
! N.B. this is the only place where d()/dA <> d()/dB
!
              G = 2.0 * D_L * F * RcosDELTA
              dGdA = 2.0 * D_L * RcosDELTA * (dFdA + 
     1                F*tanDELTA*dDELTAdA)
              dGdB = dGdA + 2.0 * F * RcosDELTA
            else
              G = 2.0 * S_L * F * RcosDELTA
              dGdB = 2.0 * S_L * RcosDELTA
     1                  *(dFdA + F * tanDELTA * dDELTAdA)
              dGdA = dGdB + 2.0 * F * RcosDELTA
            endif
          else
            G = (-1.0 + ApB * F) * RcosDELTA
            dGdA = RcosDELTA * (F - tanDELTA * dDELTAdA + ApB * dFdA
     1               + ApB * F * tanDELTA * dDELTAdA)
            dGdB = dGdA
          endif
          tmp = PsVoigt(TwoTh-DELTA+TwoTH0,TwoTH0,eta,Gamma,dPRdT
     1          ,dPRdG,dPRdE)
          sumWG = sumWG + wp(k+it) * G
          sumWRG = sumWRG + wp(k+it) * G * tmp
          sumWdGdA = sumWdGdA + wp(k+it) * dGdA
          sumWdGdB = sumWdGdB + wp(k+it) * dGdB
          sumWRdGdA = sumWRdGdA + wp(k+it) * dGdA * tmp
          sumWRdGdB = sumWRdGdB + wp(k+it) * dGdB * tmp 
          sumWGdRd2t = sumWGdRd2t + wp(k+it) * G * dPRdT
          sumWGdRdG = sumWGdRdG + wp(k+it) * G * dPRdG
          sumWGdRdE = sumWGdRdE + wp(k+it) * G * dPRdE
          sumWGdRdA = sumWGdRdA + wp(k+it) * G * dPRdT * dDELTAdA * RAD
        enddo
        if (sumWG .eq. 0.0) sumWG = 1.0
        Profval = sumWRG / sumWG
        dPRdT = sumWGdRd2t/ sumWG
        dPRdG = sumWGdRdG / sumWG
        dPRdE = sumWGdRdE / sumWG
! JPW modify here to make sum and difference be the variables:
!        dPRdS = (sumWRdGdA + sumWGdRdA) / sumWG - sumWRG *
!     1          sumWdGdA / sumWG**2
!        dPRdD = (sumWRdGdB + sumWGdRdA) / sumWG - sumWRG *
!     1          sumWdGdB / sumWG**2
        dPRdS=sumWRdGdA/sumWG
        dPRdD=sumWRdGdB/sumWG
      else   ! here for no asymmetry 
        tmp = PsVoigt(TwoTH,TwoTH0,eta,Gamma,dPRdT,dPRdG,dPRdE)
        Profval = tmp
        dPRdS = 0.0
        dPRdD = 0.0
      endif
      return
      end


      real*4 function Gauss(Pos , Pos0 , Gamma , dGdT , dGdG )

c  Return value of Gaussian at 'Pos' for peak at 'Pos0' and 'Gamma'.
c  dGdT is derivative of G wrt Pos0.
c  dGdG is derivative of G wrt Gamma.

      implicit none
      real*4 Pos , Pos0 , Gamma , dGdT , dGdG
      real*4   c / 1.6651092/
      real*4  cg / 0.939437279/
      real*4  delp , temp 

      delp = Pos - Pos0
      if ((abs(delp)*c/Gamma)**2 .gt. 40.0) then
        Gauss = 0.0
        dGdT = 0.0
        dGdG = 0.0
      else
        temp = cg * exp(-(delp * c /Gamma)**2)/Gamma
        Gauss = temp
        dGdG = temp * ( -1.0 + 2.0 * (delp * c/Gamma)**2) / Gamma
        dGdT = 2.0 * c**2 * delp * temp/Gamma**2
      endif
      return
      end

      real*4 function Lorentz(Pos , Pos0 , Gamma , dLdT , dLdG )

c  Return value of Lorentzian at 'Pos' for peak at 'Pos0' and 'Gamma'.
c  dLdT is derivative of L wrt Pos0.
c  dLdG is derivative of L wrt Gamma.
      implicit none      
      real*4 Pos , Pos0 , Gamma , dLdT , dLdG
      real*4 cl / 0.636619772/

      real*4  delp , denom

      delp = Pos - Pos0
      denom = 4.0 * delp**2 + Gamma**2
      Lorentz = cl * Gamma / denom
      dLdT = 8.0 * cl * Gamma * delp / denom**2
      dLdG = cl * (4.0 * delp**2 - Gamma**2) / denom**2
      return
      end
C                         
      real*4 function PsVoigt(TwoTH , TwoTH0 , eta , Gamma,
     1         dPRdT , dPRdG , dPRdE ) 
c
c   Returns value of Pseudo Voigt
c   Eta is the mixing coefficient between Gaussian and Lorentzian
c   Gamma is the FWHM
c   TwoTH is point at which to evaluate the profile
c   TwoTH0 is two theta value for peak
c   dPRdT is derivative of profile wrt TwoTH0
c   dPRdG is derivative of profile wrt Gamma
c   dPRdE is derivative of profile wrt Eta
      implicit none
      real*4 TwoTH , TwoTH0 , eta , Gamma
      real*4 dPRdT , dPRdG , dPRdE
      real*4  G,Gauss		! Gaussian part 
      real*4  L,Lorentz		! Lorentzian part 
      real*4 dGdT , dGdG , dLdT , dLdG
      G = Gauss(TwoTH , TwoTH0 , Gamma , dGdT , dGdG )
      L = Lorentz(TwoTH , TwoTH0 , Gamma , dLdT , dLdG )
      PsVoigt = Eta * L + (1.0 - Eta) * G
      dPRdT = Eta * dLdT + (1.0 - Eta) * dGdT
      dPRdG = Eta * dLdG + (1.0 - Eta) * dGdG
      dPRdE = L - G
      return
      end
@}




\section{Time of flight neutron}

\section{Peak widths as a function of Q}

\section{Peak fitting}

Standalone peak fitting is a handy thing to be able to do prior to indexing
a pattern, or for strain scanning type experiments.
Assuming we have a dataset containing x-values in memory and some y-values
and error bars we can attempt a peak fit. 
Initial values for the parameters are a great help.

\subsection{Fitting a single peak}

Estimation and fitting of a single peak via a single python
function. It will take a x,y,e arrays as arguements and optional
parameters for the initial estimates. Need to estimate
the initial values if they are not present.

@o profvalplusback.py
@{
@< pycopyright @>
from Numeric import *
import model
import profval
class profvalplusback(model.model):
   def __init__(self,**kwds):
      """
      Fits a single peak to some data
      Optional arguments will be estimated when computing if not present
      You are expected to "know" that keyword args are used to supply
      initial estimated values - and get those keywords from the variable
      name list
      """
      self.variables=['back','area','center','width','eta','s_l','d_l']
      self.use_asym=0
      if "s_l" in kwds.keys(): self.use_asym=1
      if "d_l" in kwds.keys(): self.use_asym=1
      model.model.__init__(self,**kwds)
      # Gradients are ycalc[0]/area for area
      #               1             for back
      self.gl={} # gradient location lookups
      self.gl['center']=1  #               ycalc[1] for center
      self.gl['width'] =2  #               ycalc[2] for width
      self.gl['eta']   =3  #               ycalc[3] for eta
      self.gl['s_l']   =4  #               ycalc[4] for s_l
      self.gl['d_l']   =5  #               ycalc[5] for d_l

   def get_variables(self):
      """
      Override the default behaviour by turning off s_l/d_l if we want to
      """
      if self.use_asym==1: return self.variables
      else:                return self.variables[0:5]

   def estimate(self,data):
      xymax=argmax(data.y)
      xymin=argmin(data.y)
      if not self.set['back']:
         self.vv[self.vd['back']]=data.y[xymin]
         self.set['back']=True
      if not self.set['center']:
         self.vv[self.vd['center']]=data.x[xymax]
         self.set['center']=True
      if not self.set['width']:
      # have to walk array
         halfway=0.5*(data.y[xymin]+data.y[xymax])
         i=xymax
         xl=data.x[0]
         while i > 0: 
            i=i-1
            if data.y[i]<halfway:
               xl=data.x[i]
               break
         i=xymax
         xh=data.x[-1]
         while i < data.x.shape[0]: 
            i=i+1
            if data.y[i]<halfway:
               xh=data.x[i]
               break
      self.vv[self.vd['width']] = abs(xl-xh)
         self.set['width']=True
      if not self.set['eta']:
         self.vv[self.vd['eta']]=0.5
         self.set['eta']=True
      if not self.set['area']:
         self.vv[self.vd['area']]=self.vv[self.vd['width']]*(data.y[xymax]-data.y[xymin])
         self.set['area']=True
      if not self.set['s_l'] and self.use_asym==1:
         self.vv[self.vd['s_l']]=0.0005
         self.set['s_l']=True
      if not self.set['d_l'] and self.use_asym==1:
         self.vv[self.vd['d_l']]=0.0005
         self.set['d_l']=True
      if self.use_asym==0:
         self.vv[self.vd['s_l']]=self.vv[self.vd['d_l']]=0.
         self.set['s_l']=True
         self.set['d_l']=True


   def compute(self,data):
      if False in self.set.values(): self.estimate(data)
      eta    = self.vv[self.vd['eta']]
      center = self.vv[self.vd['center']]
      gamma  = self.vv[self.vd['width']]
      s_l    = self.vv[self.vd['s_l']]
      d_l    = self.vv[self.vd['d_l']]
      area   = self.vv[self.vd['area']]
      #print "Computing ",eta,gamma,s_l,d_l,center,area,self.use_asym
      self.pva=profval.profval_array(eta, gamma, s_l, d_l,data.x.astype(Float32),
               center, area, self.use_asym)
      self.ycalc=self.pva[0]+self.vv[self.vd['back']]

   def gradient(self,variable):
      """
      self.pva is a tuple holding:
       ( peak_calc, dprdt , dprdg, dprde, dprds, dprdd )
      Returns the appropriate bit of it via the dictionary self.gl (gradient location)
      """
      if variable=='area':
         return (self.ycalc/self.vv[self.vd['area']]).astype(Float32)
      if variable=='back':
         return ones(self.ycalc.shape,Float32)
      return self.pva[ self.gl[variable]  ] 
@}


@o peakfit.py
@{
@< pycopyright @>

import profvalplusback
import densemodelfit
         

if __name__=="__main__":
   import epffile, powbase, mcadata, sys
   if len(sys.argv)<3:
      print "Usage: %s filename format [low] [high]"%(sys.argv[0])
      print "Formats are [epf|mca|powbase] where epf is x,y,errorbar ascii" 
      sys.exit()
   else:
      try:
         if sys.argv[2]=="powbase":
            dat=powbase.powbase(sys.argv[1])
         if sys.argv[2]=="epf":
            dat=epffile.epffile(sys.argv[1])
         if sys.argv[2]=="mca":
            dat=mcadata.mcadata(sys.argv[1])
      except:
         print "Could not read your file %s" % (sys.argv[1])
         raise
   if len(sys.argv)>3:
      lh=map(float,sys.argv[3:5])
      lh.sort()
      dat.setrange(lh)
   model=profvalplusback.profvalplusback()
   optimiser=densemodelfit.densemodelfit(model=model,data=dat)
   optimiser.refine(ncycles=1)
   while raw_input("More cycles?") != "":
      optimiser.refine(ncycles=1)
@}


\chapter{Peak center functions}

Being able to fit peak positions to raw data still leaves the user with
the problem of converting some measurement unit into some measure
of scattering vector.
In the simplest cases this is just a question of using Bragg's law, with
some more or less complicated stuff to add in depending on the experiment.

\section{Conventional fixed wavelength powder diffraction}

To a first approximation the position of a peak is given by Bragg's law, 
with a correction for a potentially imperfect setting of the zero
angle of the diffractometer.

\[ \lambda = 2 d sin(\theta - zero/2) \]

Remember that $\theta$ is half of the angle between the incident and
diffracted beams.
Other conventions for zero are possible. This one is chosen such that
if you actually measure the position of the direct beam (trivial to do
with an attenuatator), then the zero number is the angle where the
centre of the direct beam actually is.

Normally we are more likely to want $Q$ relations:

\[ Q = frac{1}{d^2} = \frac{4 \sin^2(  \theta-zero /2)}{\lambda^2} = 
                      \frac{4 \sin^2((2\theta-zero)/2)}{\lambda^2} \]
\[ \frac{dQ}{d2\theta} = \frac{8 \sin((2\theta-zero)/2) \cos((2\theta-zero)/2) }{\lambda^2} =
 - \frac{dQ}{d2\theta} \]

Note the minus sign for the zero shift. Note also this is derivative with
respect to $2\theta$ in radians. Most people use degrees for experiments
so you get an extra factor $2\pi$ in front. 

Specimen displacement errors or transparency have a particular functional
form which should be looked up and inserted into the document about here.
FIXME.

At some point we should consider refractive index corrections. Although
small these can be significant with higher resolution data. This needs
to take into account the direction of the surfaces with respect
to the directions of incoming and outgoing beams, so for small spherical
particles ($ r << \mu$) then I have a vague hope it might average
out to about zero? (FIXME)

Diffraction onto an area detector is considerably more complicated and
is going to get a section to itself later on.

\section{Fixed angle diffraction: ToF neutron and EDX}

Once again we approach the problem via Bragg's law, but this time
the angle is fixed and the wavelength is varied:

\[ \lambda = 2 d \sin(theta) \]

\subsection{ToF Neutron}

The wavelength of a neutron is given via the de'Broglie formula:

\[ \lambda = \frac{h}{p} = \frac{h}{m_n v} = \frac{h t}{m_n L} \]

where $h$ is Planck's constant, $p$ is momentum, $m_n$ is the mass
of a neutron and $v$ is the velocity. For elastic scattering, which
is all we deal with here(!), the velocity is the total flight path, $L$,
divided by the time of flight $t$.

So that we get Bragg's law in terms of neutron time of flight as:

\[ t = \frac{2 m_n L}{h} d \sin(\theta) \]

The mass of a neutron is $m_n=1.67492716\times10^{-27} kg$ and we shall
generally be working far below relativistic velocities so that no
such correction is required. Planck's contant is $h=6.62606876\times10^{-34} Js$

To relate time of flight in microseconds ($10^{-6}s$) to $d$-spacings in \AA ngstrom
we get:

\[ t [\mu s] = 505.5568 d [\AA] L \sin(\theta) + zero \]

Here zero is an offset between whatever you called time zero and the time the 
neutron left the place you measured the distance from. Due to the 
effect of a long curved guide tube there is sometimes an effectively longer 
path length for the slower neutrons, which bounce off the walls of the tube
more. Eventually in practice you have something like:

\[ t [\mu s] = 505.5568 d[\AA] [ L \sin(\theta) ] + DIFA d[\AA]^2 + zero = 
                                      DIFC d[\AA] + DIFA d[\AA]^2 + zero \]

where $DIFC = 505.5568 L \sin(\theta)$ as defined in GSAS and fullprof and
all of the three parameters are refineable.

In practice we would want to relate the time of flight to $Q=1/d^2$ and 
refine upon the angle of the counter bank and flightpath, as well as 
zero and a fudge factor. 

Best to leave it to someone working at a spallation source for now (FIXME - please!)


\section{Energy dispersive x-ray}

Despite the authors doubts about this technique (have you ever tried
to get an MCA working?) there are people wishing to fit data collected 
using a white beam of x-rays. 
We just copy the equations from the GSAS manual:
\[ E [keV] = \frac{ECONST}{2 d \sin(\theta)} + ZERO \]
where $ECONST=hc/e=12.39842 [keV \AA]$. 
The energy is related to channel via a polynomial expansion:
\[ E [keV] = \sum_{i=0}^{3} a_i c^i \]
Coefficients are determined by the biasing voltage on the detector
and the setting on the amplifier. (And the humidity in the electronics rack
I suspect.)

\section{Area detector powder diffraction}

Massive increases in detection efficiency are possible if
you just make the detector bigger and more sensitive. 
Placing a big flat plate of some sort behind the sample
and recording images gives powder data on a very fast time
scale but leads to some interesting problems of using
the data.

The images are digitised in a computer and then there is 
a question of relating the pixel position (or array indices)
to an actual postion in space and scattering vector.
Various devices are currently in use for collecting such
data and these are frequently treated using the "fit2d" 
program developed by Hammersley~\cite{fit2d} at the ESRF.

Under certain circumstances it may be beneficial to carry
out a peak intensity extraction using the raw image directly,
instead of reducing the images into 1 dimensional profiles. 
Mainly for studies of texture effects. In such cases we
will need to go through the (agonising) process of calibrating
such images.

A two step approach is forseen - firstly the intrinsic spatial 
distortion of the detector is corrected using a grid image and
then the position and orientation of the detector with respect to 
the sample and incident beam are calibrated using a standard sample.

A very nice paper about these calibrations is Stanton et al
J. Appl. Cryst (1992) 25, 549.

\subsection{Spatial distortion calibrations}

We will assume that the experimentalist has collected an image
on the detector given by a grid of holes with a well know 
spatial separation. A good experimentalist will have also taken the
same image without the grid being present for normalisation purposes.
Usually it is sufficient to place a flourescent material in the beam
as far away as possible from the detector for collecting these two
images, which should have good counting statistics.

Then we need to identify the peak positions corresponding to 
the holes. Eventually some two dimensional peak fitting might
be nice, but just using a threshold and centre of mass 
is hopefully going to give sufficient accuracy.


As an intermediate step in sorting out detector distortions we can calculate
the distortion from a fit2d spline file.
Based on the output from fit2d, it seems that the spline function used comes
from "Curve and Surface Fitting with Splines" by Paul Dierckx. There is 
code on netlib (www.netlib.org), and a brief examination suggests that
the only routine there for creating a spline function from a set of spline 
coefficients is bispev. So we just read in the information in the fit2d file
and feed it to the fortran routine from netlib, writing out the generated
function. Fortunately scipy (and pythonesrf) have the bispev routine
compiled in and ready to use.
Whether this is meaningful remains to be seen....

\subsubsection{Fit2D spline file}

An example is reproduced here:

@o spatial.spline
@{SPATIAL DISTORTION SPLINE INTERPOLATION COEFFICIENTS

  VALID REGION
 0.0000000E+00 0.0000000E+00 0.1024000E+04 0.1024000E+04

  GRID SPACING, X-PIXEL SIZE, Y-PIXEL SIZE
 0.2500000E+04 0.1620000E+03 0.1620000E+03

  X-DISTORTION
    10     9
 0.0000000E+00 0.0000000E+00 0.0000000E+00 0.0000000E+00 0.1364633E+03
 0.4858396E+03 0.1024000E+04 0.1024000E+04 0.1024000E+04 0.1024000E+04
 0.0000000E+00 0.0000000E+00 0.0000000E+00 0.0000000E+00 0.5689562E+03
 0.1024000E+04 0.1024000E+04 0.1024000E+04 0.1024000E+04
 0.3185454E+02 0.9606068E+01-0.1357207E+02 0.5400382E+01 0.1873429E+02
 0.2382785E+02 0.2834370E+01-0.1914932E+02-0.1357470E+01 0.1166279E+02
 0.1885732E+01-0.1434234E+02-0.3023303E+02-0.1834026E+02-0.7492480E+01
-0.1275525E+01-0.1536081E+01 0.2335825E+01-0.5070856E+01-0.8359183E+01
 0.1053047E+02 0.2636535E+02 0.4205279E+02 0.2190117E+02 0.9308982E+01
-0.1018729E+02 0.9855986E+01 0.3153052E+02 0.5802972E+01-0.1066609E+02

  Y-DISTORTION
    11    12
 0.0000000E+00 0.0000000E+00 0.0000000E+00 0.0000000E+00 0.1737154E+03
 0.4799059E+03 0.7900435E+03 0.1024000E+04 0.1024000E+04 0.1024000E+04
 0.1024000E+04
 0.0000000E+00 0.0000000E+00 0.0000000E+00 0.0000000E+00 0.5864732E+01
 0.5660184E+02 0.3969908E+03 0.8592832E+03 0.1024000E+04 0.1024000E+04
 0.1024000E+04 0.1024000E+04
 0.7345081E+02 0.7354512E+02 0.6850424E+02 0.4136404E+02 0.1561430E+02
 0.3190232E+02 0.2255431E+02 0.2316377E+02 0.5745498E+02 0.6135424E+02
 0.5860899E+02 0.3115641E+02 0.9689258E+01 0.3025048E+02 0.2479893E+02
 0.2321012E+02 0.2804912E+02 0.3388798E+02 0.3200676E+02 0.6025849E+01
-0.7140792E+01 0.2787918E+02 0.2903409E+02 0.2814244E+02 0.7720649E+01
 0.9234440E+01 0.5093203E+01-0.1792829E+02-0.2254351E+02 0.2528346E+02
 0.3668598E+02 0.3322346E+02 0.1690904E+02 0.1724001E+02 0.1177163E+02
-0.1184319E+02-0.1800603E+02 0.2932653E+02 0.3607935E+02 0.3693420E+02
 0.3243142E+02 0.3640476E+02 0.3562191E+02 0.9681550E+01-0.2667189E+01
 0.3414209E+02 0.3554318E+02 0.3298592E+02 0.4862249E+02 0.5061795E+02
 0.4763352E+02 0.2109127E+02 0.5088071E+01 0.3571564E+02 0.3352092E+02
 0.3274477E+02
@}

We assume the following:
\begin{itemize}
\item The valid region refers the \emph{whole} detector area
\item The grid spacing and pixel size are irrelevant to us for now
\item The two numbers after X-Distortion are the arguments nx, ny to
 the spline routine... namely the number of knots in the x and 
 y directions
\item There will be nx numbers (5 per line) in fortran format \code{5E14.7E2}
   which are the x knot positions
\item There will be ny numbers as for the nx, starting on a fresh line
\item The spline is a cubic spline having degrees kx=3 and ky=3
\item There will be (nx-kx-1)*(ny-ky-1) numbers, formatted as the knot
  positions which give the spline co-efficients
\item The Y-DISTORTION is the same as the X-DISTORTION
\end{itemize}

With those assumptions, and a careful reading of the comments in the 
top of the bispev routine which is attached in the appendix, we can
read in the file and call the subroutine, then write out the 
generated function

\subsubsection{Using Dierckx's routines from python}

Happily the folks at scipy, www.scipy.org, have already compiled 
the routine from fitpack into a python friendly version, so we
can directly use the routine from there (wahey!).
The appropriate bit of code is found in \$PYTHON/Lib/site-packages/scipy/interpolate/fitpack.py:

Sadly those same folks have made a package which is so large
and riddled with dependencies that many linux boxes appear to 
wet their pants at the very thought of installing it, so the 
appropriate bits of code were
extracted into an appendix to this document. 

\begin{verbatim}
def bisplev(x,y,tck,dx=0,dy=0):
    """Evaluate a bivariate B-spline and its derivatives.

    Description:

      Return a rank-2 array of spline function values (or spline derivative
      values) at points given by the cross-product of the rank-1 arrays x and y.
      In special cases, return an array or just a float if either x or y or
      both are floats.

    Inputs:

      x, y -- Rank-1 arrays specifying the domain over which to evaluate the
              spline or its derivative.
      tck -- A sequence of length 5 returned by bisplrep containing the knot
             locations, the coefficients, and the degree of the spline:
             [tx, ty, c, kx, ky].
      dx, dy -- The orders of the partial derivatives in x and y respectively.

    Outputs: (vals, )

      vals -- The B-pline or its derivative evaluated over the set formed by
              the cross-product of x and y.
    """
\end{verbatim}

So to do the job of reading a fi2d spline file in python, we just need
to set up the tck tuple from what is in the file and ask for the 
x,y array we want back.

@D readfit2dspline
@{
from Numeric import *
import bisplev
   
def readfit2dfloats(fp,n,debug=0):
   """
   Interprets a 5E14.7 formatted fortran line
   """
   vals=[]
   j=0
   while j < n:
      i=0
      line=fp.readline()
      while i < 5*14:
         if(debug):print line
         if(debug):print i,line[i:i+14]
         vals.append(float(line[i:i+14]) )
         j=j+1
         i=i+14
         if j == n: break
      i=i+1
   return vals


# read the fit2d array into a tck tuple
def readfit2dspline(name):
   """
   Reads a fit2d spline file into a scipy/fitpack tuple, tck
   A fairly long and dull routine...
   """
   kx=3
   ky=3
   fp=open(name,"r")
   line=fp.readline() # SPATIAL DISTORTION SPLINE INTERPOLATION COEFFICIENTS
   if line[:7] != "SPATIAL":
      raise SyntaxError, name+": file does not seem to be a fit2d spline file"
   line=fp.readline() # BLANK LINE
   line=fp.readline() # VALID REGION
   line=fp.readline() # the actual valid region, assume xmin,ymin,xmax,ymax
   vals=line.split()
   xmin=float(vals[0])
   ymin=float(vals[1])
   xmax=float(vals[3])
   ymax=float(vals[3])
   line=fp.readline() # BLANK
   line=fp.readline() # GRID SPACING, X-PIXEL SIZE, Y-PIXEL SIZE
   line=fp.readline()
   vals=line.split()
   gs=float(vals[0])
   xps=float(vals[1])
   yps=float(vals[2])
   line=fp.readline() # BLANK
   line=fp.readline() # X-DISTORTION
   line=fp.readline() # two integers nx1,ny1
   vals=line.split()
   nx1=int(vals[0])
   ny1=int(vals[1])
   # Now follow fit2d formatted line 5E14.7
   tx1=array(readfit2dfloats(fp,nx1),Float32)
   ty1=array(readfit2dfloats(fp,ny1),Float32)
   c1 =array(readfit2dfloats(fp,(nx1-4)*(ny1-4)),Float32)
   line=fp.readline() #BLANK
   line=fp.readline() # Y-DISTORTION
   line=fp.readline() # two integers nx2, ny2
   vals=line.split()
   nx2=int(vals[0])
   ny2=int(vals[1])
   tx2=array(readfit2dfloats(fp,nx2),Float32)
   ty2=array(readfit2dfloats(fp,ny2),Float32)
   c2 =array(readfit2dfloats(fp,(nx2-4)*(ny2-4)),Float32)
   fp.close()
   tck1=(tx1,ty1,c1,kx,ky)
   tck2=(tx2,ty2,c2,kx,ky)
   return tck1,tck2,xmin,xmax,ymin,ymax
@}

@D splinetoimages
@{@<readfit2dspline@>
def splinetoimages(name):
   tck1,tck2,xmin,xmax,ymin,ymax=readfit2dspline(name)
   x=arange(xmin,xmax,1.)
   y=arange(ymin,ymax,1.)
   im1=transpose(bisplev.bisplev(x,y,tck1)) # transpose as bisplev/fit2d speak fortran
   im2=transpose(bisplev.bisplev(x,y,tck2))
   return im1, im2

def splinetopixelcoords(name):
   im1,im2=splinetoimages(name)
   iarray=fromfunction(lambda i,j: j , im1.shape )
   jarray=fromfunction(lambda i,j: i , im2.shape )
   return im1+iarray, im2+jarray

class fit2dsplinefunc:
   def __init__(self,filename=None):
      if filename!=None:
         self.tck1,self.tck2,self.xmin,self.xmax,\
                               self.ymin,self.ymax=readfit2dspline(filename)
         self.getdistortionimages()
      else: print "The spline object exists but is empty"
   def getxy(self,xy):
      """ 
      xy is an array of peak positions with dims npeaks*2
      Due to the inferface to bisplev we compute npeaks**2 and take the diagonal
      of the resulting matrix. Lovely.
      """
      xnew = diagonal(bisplev.bisplev(xy[:,0],xy[:,1],tck1))
      ynew = diagonal(bisplev.bisplev(xy[:,0],xy[:,1],tck2))
      return transpose(array([xnew,ynew],Float))
   def getdistortionimages(self):
      """
      Returns the distortion at each pixel position in the image
      """
      x=arange(self.xmin,self.xmax,1.)
      y=arange(self.ymin,self.ymax,1.)
      self.xd=transpose(bisplev.bisplev(x,y,self.tck1)) 
            # transpose as bisplev/fit2d speak fortran
      self.yd=transpose(bisplev.bisplev(x,y,self.tck2))
   def getpixelcoordinateimages(self):
      """
      Returns the floating point numbers which could replace i,j 
      as pixel positions when looking at images
      """
      if self.xd==None: self.getdistortionimages()
      return self.xd+fromfunction(lambda i,j : j, self.xd.shape), \
             self.yd+fromfunction(lambda i,j : i, self.xd.shape)

   def unwarp(self,image):
      """
      Use the xd and yd distortion images from the routine splinetoimages (reading a 
      fit2d spline file to unwarp an image). 
      Assumes the image has come from reading a bytestream so that it is laid out in
      memory in the "normal" way.
      Current version is for debugging only - it just takes the floor of the unwarped
      pixel value for deciding which bin to go into. 
      Proper rebinning should be done at some later date. 
      """
      newim=zeros(image.shape,Float32)  # To hold the result
      ind=array(range(multiply.reduce(image.shape)),Int32) # For indexing the flat array
      # ind[x,y]=ind.flat[x*ind.shape[1]+y]
      # x=i+yd[i,j] to index image[x,y] 
      # Use ravel instead of flat in case array is not contiguous
      ind=ind+ravel(floor(self.yd).astype(Int32))*image.shape[1]
      ind=ind+ravel(floor(self.xd).astype(Int32))
      ind=where(ind<0,0,ind)
      ind=where(ind>=ind.shape[0],ind.shape[0]-1,ind)
      #   print maximum.reduce(ind),minimum.reduce(ind)
      put(newim,ind,image) # This just take one pixel (the floor)
                        # it ought to also add in the ceil value and do proper binning
      return newim

@}

That covers most of the functionality of spline files that we wish to persue. There is
also a problem of testing, which means we would like to check that pixel positions
are coming out about right. By doing an "unwarp" we can check if they do or not, 
bearing in mind that I don't want to get into the messy details of splitting pixels
correctly.

@D unwarp
@{
def unwarp(image,xd,yd):
   """
   Use the xd and yd distortion images from the routine splinetoimages (reading a 
   fit2d spline file to unwarp an image). 
   Assumes the image has come from reading a bytestream so that it is laid out in
   memory in the "normal" way.
   Current version is for debugging only - it just takes the floor of the unwarped
   pixel value for deciding which bin to go into. 
   Proper rebinning should be done at some later date. 
   """
   newim=zeros(image.shape,Float32)  # To hold the result
   ind=array(range(multiply.reduce(image.shape)),Int32) # For indexing the flat array
   # ind[x,y]=ind.flat[x*ind.shape[1]+y]
   # x=i+yd[i,j] to index image[x,y] 
   # Use ravel instead of flat in case array is not contiguous
   ind=ind+ravel(floor(yd).astype(Int32))*image.shape[1]
   ind=ind+ravel(floor(xd).astype(Int32))
   ind=where(ind<0,0,ind)
   ind=where(ind>=ind.shape[0],ind.shape[0]-1,ind)
#   print maximum.reduce(ind),minimum.reduce(ind)
   put(newim,ind,image) # This just take one pixel (the floor)
                        # it ought to also add in the ceil value and do proper binning
   return newim
@}


Code to test the unwarping:

@O unwarp.py
@{
import imagereaders 
@<splinetoimages@>
@<unwarp@>
if __name__=="__main__":
   import sys
   from Numeric import *
   if len(sys.argv)<4:
      print "Usage: ",sys.argv[0]," imagefile splinefile outfile"
      sys.exit()
   infile=sys.argv[1]
   splinefile=sys.argv[2]
   outfile=sys.argv[3]
   xd,yd=splinetoimages(splinefile)
   if infile[-3:]=="edf":
      h,inarray=readid11edf(infile)
   else:
      h,inarray=readbruker(infile)
   outarray=unwarp(inarray,xd,yd)
   f=open(outfile,"wb")
   f.write(outarray.tostring())
   f.close()
   print "Wrote ",outfile," with dimensions ",outarray.shape,\
         " and type ",outarray.typecode()
   
@}

\subsubsection{bisplev 2 d spline source}

\subsection{Spline routines from scipy/fitpack}

Due to the overwhelming difficulties with getting the entire scipy package to build
and install under linux, or downloading some massive python installation for windows,
the useful bits for doing spline functions have been stripped out and reproduced
here. They say in their FAQ:
\begin{verbatim}
The intent of SciPy's licensing is that it is free for both commercial 
and non-commercial use, just don't sue us. For now we're distributing it
under the BSD license. We're not open source licensing wizards (OSLW), and 
its final format is not yet known. We'd like to use the Python license, 
but understand that it isn't really written to protect a corporation (which we are).
\end{verbatim}
Which seems fair enough and very generous of them. They also say it is "Open Source". 

Appears that I only need to include their license:

@o bisplev_license.txt
@{
This license applies to the interface between python and the spline fitting
functions from Paul Dierckx's fortran fitpack routines.

Copyright (c) 2001, 2002 Enthought, Inc.

All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

  a. Redistributions of source code must retain the above copyright notice,
     this list of conditions and the following disclaimer.
  b. Redistributions in binary form must reproduce the above copyright
     notice, this list of conditions and the following disclaimer in the
     documentation and/or other materials provided with the distribution.
  c. Neither the name of the Enthought nor the names of its contributors
     may be used to endorse or promote products derived from this software
     without specific prior written permission.


THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE REGENTS OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
DAMAGE.
@}


@O bisplev.py
@{      

# THIS CODE COMES FROM THE SCIPY PACKAGE AT www.scipy.org
#
# IT HAS BEEN COPIED TO SAVE YOU HAVING TO INSTALL ALL THE 
# OTHER STUFF FROM THERE. 
#
import _splines
from Numeric import *

def myasarray(a):
    if type(a) in [type(1.0),type(1L),type(1),type(1j)]:
        return asarray([a])
    elif type(a) is ArrayType and len(a)==1:
        # Takes care of mapping array(number) to array([number])
        return asarray([a[0]])
    else:
        return asarray(a)

def bisplev(x,y,tck,dx=0,dy=0):
   """Evaluate a bivariate B-spline and its derivatives.
   Description:
     Return a rank-2 array of spline function values (or spline derivative
     values) at points given by the cross-product of the rank-1 arrays x and y.
     In special cases, return an array or just a float if either x or y or
     both are floats.
   Inputs:
     x, y -- Rank-1 arrays specifying the domain over which to evaluate the
             spline or its derivative.
     tck -- A sequence of length 5 returned by bisplrep containing the knot
            locations, the coefficients, and the degree of the spline:
            [tx, ty, c, kx, ky].
     dx, dy -- The orders of the partial derivatives in x and y respectively.
   Outputs: (vals, )
     vals -- The B-pline or its derivative evaluated over the set formed by
             the cross-product of x and y.
   """
   tx,ty,c,kx,ky=tck
   if not (0<=dx<kx): raise ValueError,"0<=dx=%d<kx=%d must hold"%(dx,kx)
   if not (0<=dy<ky): raise ValueError,"0<=dy=%d<ky=%d must hold"%(dy,ky)
   x,y=map(myasarray,[x,y])
   if (len(x.shape) != 1) or (len(y.shape) != 1):
       raise ValueError, "First two entries should be rank-1 arrays."
   z,ier=_splines._bispev(tx,ty,c,kx,ky,x,y,dx,dy)
   if ier==10: raise ValueError,"Invalid input data"
   if ier: raise TypeError,"An error occurred"
   z.shape=len(x),len(y)
   if len(z)>1: return z
   if len(z[0])>1: return z[0]
   return z[0][0]



_surfit_cache = {'tx': array([],'d'),'ty': array([],'d'),
                 'wrk': array([],'d'), 'iwrk':array([],'i')}

def bisplrep(x,y,z,w=None,xb=None,xe=None,yb=None,ye=None,kx=3,ky=3,task=0,s=None,
             eps=1e-16,tx=None,ty=None,full_output=0,nxest=None,nyest=None,quiet=1):
    """Find a bivariate B-spline representation of a surface.

    Description:

      Given a set of data points (x[i], y[i], z[i]) representing a surface
      z=f(x,y), compute a B-spline representation of the surface.

    Inputs:

      x, y, z -- Rank-1 arrays of data points.
      w -- Rank-1 array of weights. By default w=ones(len(x)).
      xb, xe -- End points of approximation interval in x. 
      yb, ye -- End points of approximation interval in y.
                By default xb, xe, yb, ye = x[0], x[-1], y[0], y[-1]
      kx, ky -- The degrees of the spline (1 <= kx, ky <= 5).  Third order
                (kx=ky=3) is recommended.
      task -- If task=0, find knots in x and y and coefficients for a given
                smoothing factor, s.
              If task=1, find knots and coefficients for another value of the
                smoothing factor, s.  bisplrep must have been previously called
                with task=0 or task=1.
              If task=-1, find coefficients for a given set of knots tx, ty.
      s -- A non-negative smoothing factor.  If weights correspond to the inverse
           of the standard-deviation of the errors in z, then a good s-value
           should be found in the range (m-sqrt(2*m),m+sqrt(2*m)) where m=len(x)
      eps -- A threshold for determining the effective rank of an over-determined
             linear system of equations (0 < eps < 1) --- not likely to need
             changing.
      tx, ty -- Rank-1 arrays of the knots of the spline for task=-1
      full_output -- Non-zero to return optional outputs.
      nxest, nyest -- Over-estimates of the total number of knots.  If None then
                      nxest = max(kx+sqrt(m/2),2*kx+3),
                      nyest = max(ky+sqrt(m/2),2*ky+3)
      quiet -- Non-zero to suppress printing of messages.

    Outputs: (tck, {fp, ier, msg})

      tck -- A list [tx, ty, c, kx, ky] containing the knots (tx, ty) and
             coefficients (c) of the bivariate B-spline representation of the
             surface along with the degree of the spline.

      fp -- The weighted sum of squared residuals of the spline approximation.
      ier -- An integer flag about splrep success.  Success is indicated if
             ier<=0. If ier in [1,2,3] an error occurred but was not raised.
             Otherwise an error is raised.
      msg -- A message corresponding to the integer flag, ier.

    Remarks:

      SEE bisplev to evaluate the value of the B-spline given its tck
      representation.           
    """
    x,y,z=map(myasarray,[x,y,z])
    x,y,z=map(ravel,[x,y,z])  # ensure 1-d arrays.
    m=len(x)
    if not (m==len(y)==len(z)): raise TypeError, 'len(x)==len(y)==len(z) must hold.'  
    if w is None: w=ones(m,'d')
    else: w=myasarray(w)
    if not len(w) == m: raise TypeError,' len(w)=%d is not equal to m=%d'%(len(w),m)
    if xb is None: xb=x[0]
    if xe is None: xe=x[-1]
    if yb is None: yb=y[0]
    if ye is None: ye=y[-1]
    if not (-1<=task<=1): raise TypeError, 'task must be either -1,0, or 1'
    if s is None: s=m-sqrt(2*m)
    if tx is None and task==-1: raise TypeError, 'Knots_x must be given for task=-1'
    if tx is not None: _curfit_cache['tx']=myasarray(tx)
    nx=len(_surfit_cache['tx'])
    if ty is None and task==-1: raise TypeError, 'Knots_y must be given for task=-1'
    if ty is not None: _curfit_cache['ty']=myasarray(ty)
    ny=len(_surfit_cache['ty'])
    if task==-1 and nx<2*kx+2:
        raise TypeError, 'There must be at least 2*kx+2 knots_x for task=-1'
    if task==-1 and ny<2*ky+2:
        raise TypeError, 'There must be at least 2*ky+2 knots_x for task=-1'
    if not ((1<=kx<=5) and (1<=ky<=5)): 
        raise TypeError, \
       'Given degree of the spline (kx,ky=%d,%d) is not supported. (1<=k<=5)'%(kx,ky)
    if m<(kx+1)*(ky+1): raise TypeError, 'm>=(kx+1)(ky+1) must hold'
    if nxest is None: nxest=kx+sqrt(m/2)
    if nyest is None: nyest=ky+sqrt(m/2)
    nxest,nyest=max(nxest,2*kx+3),max(nyest,2*ky+3)
    if task>=0 and s==0:
        nxest=int(kx+sqrt(3*m))
        nyest=int(ky+sqrt(3*m))
    if task==-1:
        _surfit_cache['tx']=myasarray(tx)
        _surfit_cache['ty']=myasarray(ty)
    tx,ty=_surfit_cache['tx'],_surfit_cache['ty']
    wrk=_surfit_cache['wrk']
    iwrk=_surfit_cache['iwrk']
    u,v,km,ne=nxest-kx-1,nyest-ky-1,max(kx,ky)+1,max(nxest,nyest)
    bx,by=kx*v+ky+1,ky*u+kx+1
    b1,b2=bx,bx+v-ky
    if bx>by: b1,b2=by,by+u-kx
    lwrk1=u*v*(2+b1+b2)+2*(u+v+km*(m+ne)+ne-kx-ky)+b2+1
    lwrk2=u*v*(b2+1)+b2
    tx,ty,c,o = _fitpack._surfit(x,y,z,w,xb,xe,yb,ye,kx,ky,task,s,eps,
                                   tx,ty,nxest,nyest,wrk,lwrk1,lwrk2)
    _curfit_cache['tx']=tx
    _curfit_cache['ty']=ty
    _curfit_cache['wrk']=o['wrk']
    ier,fp=o['ier'],o['fp']
    tck=[tx,ty,c,kx,ky]
    if ier<=0 and not quiet:
        print _iermess2[ier][0]
        print "\tkx,ky=%d,%d nx,ny=%d,%d m=%d fp=%f s=%f"%(kx,ky,len(tx),
                                                           len(ty),m,fp,s)
    ierm=min(11,max(-3,ier))
    if ierm>0 and not full_output:
        if ier in [1,2,3,4,5]:
            print "Warning: "+_iermess2[ierm][0]
            print "\tkx,ky=%d,%d nx,ny=%d,%d m=%d fp=%f s=%f"%(kx,ky,len(tx),
                                                           len(ty),m,fp,s)
        else:
            try:
                raise _iermess2[ierm][1],_iermess2[ierm][0]
            except KeyError:
                raise _iermess2['unknown'][1],_iermess2['unknown'][0]
    if full_output:
        try:
            return tck,fp,ier,_iermess2[ierm][0]
        except KeyError:
            return tck,fp,ier,_iermess2['unknown'][0]
    else:
        return tck

@}


@o bispev.f
@{
      subroutine bispev(tx,nx,ty,ny,c,kx,ky,x,mx,y,my,z,wrk,lwrk,
     * iwrk,kwrk,ier)
c  subroutine bispev evaluates on a grid (x(i),y(j)),i=1,...,mx; j=1,...
c  ,my a bivariate spline s(x,y) of degrees kx and ky, given in the
c  b-spline representation.
c
c  calling sequence:
c     call bispev(tx,nx,ty,ny,c,kx,ky,x,mx,y,my,z,wrk,lwrk,
c    * iwrk,kwrk,ier)
c
c  input parameters:
c   tx    : real array, length nx, which contains the position of the
c           knots in the x-direction.
c   nx    : integer, giving the total number of knots in the x-direction
c   ty    : real array, length ny, which contains the position of the
c           knots in the y-direction.
c   ny    : integer, giving the total number of knots in the y-direction
c   c     : real array, length (nx-kx-1)*(ny-ky-1), which contains the
c           b-spline coefficients.
c   kx,ky : integer values, giving the degrees of the spline.
c   x     : real array of dimension (mx).
c           before entry x(i) must be set to the x co-ordinate of the
c           i-th grid point along the x-axis.
c           tx(kx+1)<=x(i-1)<=x(i)<=tx(nx-kx), i=2,...,mx.
c   mx    : on entry mx must specify the number of grid points along
c           the x-axis. mx >=1.
c   y     : real array of dimension (my).
c           before entry y(j) must be set to the y co-ordinate of the
c           j-th grid point along the y-axis.
c           ty(ky+1)<=y(j-1)<=y(j)<=ty(ny-ky), j=2,...,my.
c   my    : on entry my must specify the number of grid points along
c           the y-axis. my >=1.
c   wrk   : real array of dimension lwrk. used as workspace.
c   lwrk  : integer, specifying the dimension of wrk.
c           lwrk >= mx*(kx+1)+my*(ky+1)
c   iwrk  : integer array of dimension kwrk. used as workspace.
c   kwrk  : integer, specifying the dimension of iwrk. kwrk >= mx+my.
c
c  output parameters:
c   z     : real array of dimension (mx*my).
c           on succesful exit z(my*(i-1)+j) contains the value of s(x,y)
c           at the point (x(i),y(j)),i=1,...,mx;j=1,...,my.
c   ier   : integer error flag
c    ier=0 : normal return
c    ier=10: invalid input data (see restrictions)
c
c  restrictions:
c   mx >=1, my >=1, lwrk>=mx*(kx+1)+my*(ky+1), kwrk>=mx+my
c   tx(kx+1) <= x(i-1) <= x(i) <= tx(nx-kx), i=2,...,mx
c   ty(ky+1) <= y(j-1) <= y(j) <= ty(ny-ky), j=2,...,my
c
c  other subroutines required:
c    fpbisp,fpbspl
c
c  references :
c    de boor c : on calculating with b-splines, j. approximation theory
c                6 (1972) 50-62.
c    cox m.g.  : the numerical evaluation of b-splines, j. inst. maths
c                applics 10 (1972) 134-149.
c    dierckx p. : curve and surface fitting with splines, monographs on
c                 numerical analysis, oxford university press, 1993.
c
c  author :
c    p.dierckx
c    dept. computer science, k.u.leuven
c    celestijnenlaan 200a, b-3001 heverlee, belgium.
c    e-mail : Paul.Dierckx@@cs.kuleuven.ac.be
c
c  latest update : march 1987
c
c  ..scalar arguments..
      integer nx,ny,kx,ky,mx,my,lwrk,kwrk,ier
c  ..array arguments..
      integer iwrk(kwrk)
      real*8 tx(nx),ty(ny),c((nx-kx-1)*(ny-ky-1)),x(mx),y(my),z(mx*my),
     * wrk(lwrk)
c  ..local scalars..
      integer i,iw,lwest
c  ..
c  before starting computations a data check is made. if the input data
c  are invalid control is immediately repassed to the calling program.
      ier = 10
      lwest = (kx+1)*mx+(ky+1)*my
      if(lwrk.lt.lwest) go to 100
      if(kwrk.lt.(mx+my)) go to 100
      if(mx-1) 100,30,10
  10  do 20 i=2,mx
        if(x(i).lt.x(i-1)) go to 100
  20  continue
  30  if(my-1) 100,60,40
  40  do 50 i=2,my
        if(y(i).lt.y(i-1)) go to 100
  50  continue
  60  ier = 0
      iw = mx*(kx+1)+1
      call fpbisp(tx,nx,ty,ny,c,kx,ky,x,mx,y,my,z,wrk(1),wrk(iw),
     * iwrk(1),iwrk(mx+1))
 100  return
      end
      subroutine fpbisp(tx,nx,ty,ny,c,kx,ky,x,mx,y,my,z,wx,wy,lx,ly)
c  ..scalar arguments..
      integer nx,ny,kx,ky,mx,my
c  ..array arguments..
      integer lx(mx),ly(my)
      real*8 tx(nx),ty(ny),c((nx-kx-1)*(ny-ky-1)),x(mx),y(my),z(mx*my),
     * wx(mx,kx+1),wy(my,ky+1)
c  ..local scalars..
      integer kx1,ky1,l,l1,l2,m,nkx1,nky1
      real*8 arg,sp,tb,te
c  ..local arrays..
      real*8 h(6)
c  ..subroutine references..
c    fpbspl
c  ..
      kx1 = kx+1
      nkx1 = nx-kx1
      tb = tx(kx1)
      te = tx(nkx1+1)
      l = kx1
      l1 = l+1
      do 40 i=1,mx
        arg = x(i)
        if(arg.lt.tb) arg = tb
        if(arg.gt.te) arg = te
  10    if(arg.lt.tx(l1) .or. l.eq.nkx1) go to 20
        l = l1
        l1 = l+1
        go to 10
  20    call fpbspl(tx,nx,kx,arg,l,h)
        lx(i) = l-kx1
        do 30 j=1,kx1
          wx(i,j) = h(j)
  30    continue
  40  continue
      ky1 = ky+1
      nky1 = ny-ky1
      tb = ty(ky1)
      te = ty(nky1+1)
      l = ky1
      l1 = l+1
      do 80 i=1,my
        arg = y(i)
        if(arg.lt.tb) arg = tb
        if(arg.gt.te) arg = te
  50    if(arg.lt.ty(l1) .or. l.eq.nky1) go to 60
        l = l1
        l1 = l+1
        go to 50
  60    call fpbspl(ty,ny,ky,arg,l,h)
        ly(i) = l-ky1
        do 70 j=1,ky1
          wy(i,j) = h(j)
  70    continue
  80  continue
      m = 0
      do 130 i=1,mx
        l = lx(i)*nky1
        do 90 i1=1,kx1
          h(i1) = wx(i,i1)
  90    continue
        do 120 j=1,my
          l1 = l+ly(j)
          sp = 0.
          do 110 i1=1,kx1
            l2 = l1
            do 100 j1=1,ky1
              l2 = l2+1
              sp = sp+c(l2)*h(i1)*wy(j,j1)
 100        continue
            l1 = l1+nky1
 110      continue
          m = m+1
          z(m) = sp
 120    continue
 130  continue
      return
      end

      subroutine fpbspl(t,n,k,x,l,h)
c  subroutine fpbspl evaluates the (k+1) non-zero b-splines of
c  degree k at t(l) <= x < t(l+1) using the stable recurrence
c  relation of de boor and cox.
c  ..
c  ..scalar arguments..
      real*8 x
      integer n,k,l
c  ..array arguments..
      real*8 t(n),h(6)
c  ..local scalars..
      real*8 f,one
      integer i,j,li,lj
c  ..local arrays..
      real*8 hh(5)
c  ..
      one = 0.1d+01
      h(1) = one
      do 20 j=1,k
        do 10 i=1,j
          hh(i) = h(i)
  10    continue
        h(1) = 0.0d0
        do 20 i=1,j
          li = l+i
          lj = li-j
          f = hh(i)/(t(li)-t(lj))
          h(i) = h(i)+f*(t(li)-x)
          h(i+1) = f*(x-t(lj))
  20  continue
      return
      end
      subroutine parder(tx,nx,ty,ny,c,kx,ky,nux,nuy,x,mx,y,my,z,
     * wrk,lwrk,iwrk,kwrk,ier)
c  subroutine parder evaluates on a grid (x(i),y(j)),i=1,...,mx; j=1,...
c  ,my the partial derivative ( order nux,nuy) of a bivariate spline
c  s(x,y) of degrees kx and ky, given in the b-spline representation.
c
c  calling sequence:
c     call parder(tx,nx,ty,ny,c,kx,ky,nux,nuy,x,mx,y,my,z,wrk,lwrk,
c    * iwrk,kwrk,ier)
c
c  input parameters:
c   tx    : real array, length nx, which contains the position of the
c           knots in the x-direction.
c   nx    : integer, giving the total number of knots in the x-direction
c   ty    : real array, length ny, which contains the position of the
c           knots in the y-direction.
c   ny    : integer, giving the total number of knots in the y-direction
c   c     : real array, length (nx-kx-1)*(ny-ky-1), which contains the
c           b-spline coefficients.
c   kx,ky : integer values, giving the degrees of the spline.
c   nux   : integer values, specifying the order of the partial
c   nuy     derivative. 0<=nux<kx, 0<=nuy<ky.
c   x     : real array of dimension (mx).
c           before entry x(i) must be set to the x co-ordinate of the
c           i-th grid point along the x-axis.
c           tx(kx+1)<=x(i-1)<=x(i)<=tx(nx-kx), i=2,...,mx.
c   mx    : on entry mx must specify the number of grid points along
c           the x-axis. mx >=1.
c   y     : real array of dimension (my).
c           before entry y(j) must be set to the y co-ordinate of the
c           j-th grid point along the y-axis.
c           ty(ky+1)<=y(j-1)<=y(j)<=ty(ny-ky), j=2,...,my.
c   my    : on entry my must specify the number of grid points along
c           the y-axis. my >=1.
c   wrk   : real array of dimension lwrk. used as workspace.
c   lwrk  : integer, specifying the dimension of wrk.
c           lwrk >= mx*(kx+1-nux)+my*(ky+1-nuy)+(nx-kx-1)*(ny-ky-1)
c   iwrk  : integer array of dimension kwrk. used as workspace.
c   kwrk  : integer, specifying the dimension of iwrk. kwrk >= mx+my.
c
c  output parameters:
c   z     : real array of dimension (mx*my).
c           on succesful exit z(my*(i-1)+j) contains the value of the
c           specified partial derivative of s(x,y) at the point
c           (x(i),y(j)),i=1,...,mx;j=1,...,my.
c   ier   : integer error flag
c    ier=0 : normal return
c    ier=10: invalid input data (see restrictions)
c
c  restrictions:
c   mx >=1, my >=1, 0 <= nux < kx, 0 <= nuy < ky, kwrk>=mx+my
c   lwrk>=mx*(kx+1-nux)+my*(ky+1-nuy)+(nx-kx-1)*(ny-ky-1),
c   tx(kx+1) <= x(i-1) <= x(i) <= tx(nx-kx), i=2,...,mx
c   ty(ky+1) <= y(j-1) <= y(j) <= ty(ny-ky), j=2,...,my
c
c  other subroutines required:
c    fpbisp,fpbspl
c
c  references :
c    de boor c : on calculating with b-splines, j. approximation theory
c                6 (1972) 50-62.
c   dierckx p. : curve and surface fitting with splines, monographs on
c                numerical analysis, oxford university press, 1993.
c
c  author :
c    p.dierckx
c    dept. computer science, k.u.leuven
c    celestijnenlaan 200a, b-3001 heverlee, belgium.
c    e-mail : Paul.Dierckx@@cs.kuleuven.ac.be
c
c  latest update : march 1989
c
c  ..scalar arguments..
      integer nx,ny,kx,ky,nux,nuy,mx,my,lwrk,kwrk,ier
c  ..array arguments..
      integer iwrk(kwrk)
      real*8 tx(nx),ty(ny),c((nx-kx-1)*(ny-ky-1)),x(mx),y(my),z(mx*my),
     * wrk(lwrk)
c  ..local scalars..
      integer i,iwx,iwy,j,kkx,kky,kx1,ky1,lx,ly,lwest,l1,l2,m,m0,m1,
     * nc,nkx1,nky1,nxx,nyy
      real*8 ak,fac
c  ..
c  before starting computations a data check is made. if the input data
c  are invalid control is immediately repassed to the calling program.
      ier = 10
      kx1 = kx+1
      ky1 = ky+1
      nkx1 = nx-kx1
      nky1 = ny-ky1
      nc = nkx1*nky1
      if(nux.lt.0 .or. nux.ge.kx) go to 400
      if(nuy.lt.0 .or. nuy.ge.ky) go to 400
      lwest = nc +(kx1-nux)*mx+(ky1-nuy)*my
      if(lwrk.lt.lwest) go to 400
      if(kwrk.lt.(mx+my)) go to 400
      if(mx-1) 400,30,10
  10  do 20 i=2,mx
        if(x(i).lt.x(i-1)) go to 400
  20  continue
  30  if(my-1) 400,60,40
  40  do 50 i=2,my
        if(y(i).lt.y(i-1)) go to 400
  50  continue
  60  ier = 0
      nxx = nkx1
      nyy = nky1
      kkx = kx
      kky = ky
c  the partial derivative of order (nux,nuy) of a bivariate spline of
c  degrees kx,ky is a bivariate spline of degrees kx-nux,ky-nuy.
c  we calculate the b-spline coefficients of this spline
      do 70 i=1,nc
        wrk(i) = c(i)
  70  continue
      if(nux.eq.0) go to 200
      lx = 1
      do 100 j=1,nux
        ak = kkx
        nxx = nxx-1
        l1 = lx
        m0 = 1
        do 90 i=1,nxx
          l1 = l1+1
          l2 = l1+kkx
          fac = tx(l2)-tx(l1)
          if(fac.le.0.) go to 90
          do 80 m=1,nyy
            m1 = m0+nyy
            wrk(m0) = (wrk(m1)-wrk(m0))*ak/fac
            m0  = m0+1
  80      continue
  90    continue
        lx = lx+1
        kkx = kkx-1
 100  continue
 200  if(nuy.eq.0) go to 300
      ly = 1
      do 230 j=1,nuy
        ak = kky
        nyy = nyy-1
        l1 = ly
        do 220 i=1,nyy
          l1 = l1+1
          l2 = l1+kky
          fac = ty(l2)-ty(l1)
          if(fac.le.0.) go to 220
          m0 = i
          do 210 m=1,nxx
            m1 = m0+1
            wrk(m0) = (wrk(m1)-wrk(m0))*ak/fac
            m0  = m0+nky1
 210      continue
 220    continue
        ly = ly+1
        kky = kky-1
 230  continue
      m0 = nyy
      m1 = nky1
      do 250 m=2,nxx
        do 240 i=1,nyy
          m0 = m0+1
          m1 = m1+1
          wrk(m0) = wrk(m1)
 240    continue
        m1 = m1+nuy
 250  continue
c  we partition the working space and evaluate the partial derivative
 300  iwx = 1+nxx*nyy
      iwy = iwx+mx*(kx1-nux)
      call fpbisp(tx(nux+1),nx-2*nux,ty(nuy+1),ny-2*nuy,wrk,kkx,kky,
     * x,mx,y,my,z,wrk(iwx),wrk(iwy),iwrk(1),iwrk(mx+1))
 400  return
      end
c
      subroutine surfit(iopt,m,x,y,z,w,xb,xe,yb,ye,kx,ky,s,nxest,nyest,
     *  nmax,eps,nx,tx,ny,ty,c,fp,wrk1,lwrk1,wrk2,lwrk2,iwrk,kwrk,ier)
c given the set of data points (x(i),y(i),z(i)) and the set of positive
c numbers w(i),i=1,...,m, subroutine surfit determines a smooth bivar-
c iate spline approximation s(x,y) of degrees kx and ky on the rect-
c angle xb <= x <= xe, yb <= y <= ye.
c if iopt = -1 surfit calculates the weighted least-squares spline
c according to a given set of knots.
c if iopt >= 0 the total numbers nx and ny of these knots and their
c position tx(j),j=1,...,nx and ty(j),j=1,...,ny are chosen automatic-
c ally by the routine. the smoothness of s(x,y) is then achieved by
c minimalizing the discontinuity jumps in the derivatives of s(x,y)
c across the boundaries of the subpanels (tx(i),tx(i+1))*(ty(j),ty(j+1).
c the amounth of smoothness is determined by the condition that f(p) =
c sum ((w(i)*(z(i)-s(x(i),y(i))))**2) be <= s, with s a given non-neg-
c ative constant, called the smoothing factor.
c the fit is given in the b-spline representation (b-spline coefficients
c c((ny-ky-1)*(i-1)+j),i=1,...,nx-kx-1;j=1,...,ny-ky-1) and can be eval-
c uated by means of subroutine bispev.
c
c calling sequence:
c     call surfit(iopt,m,x,y,z,w,xb,xe,yb,ye,kx,ky,s,nxest,nyest,
c    *  nmax,eps,nx,tx,ny,ty,c,fp,wrk1,lwrk1,wrk2,lwrk2,iwrk,kwrk,ier)
c
c parameters:
c  iopt  : integer flag. on entry iopt must specify whether a weighted
c          least-squares spline (iopt=-1) or a smoothing spline (iopt=0
c          or 1) must be determined.
c          if iopt=0 the routine will start with an initial set of knots
c          tx(i)=xb,tx(i+kx+1)=xe,i=1,...,kx+1;ty(i)=yb,ty(i+ky+1)=ye,i=
c          1,...,ky+1. if iopt=1 the routine will continue with the set
c          of knots found at the last call of the routine.
c          attention: a call with iopt=1 must always be immediately pre-
c                     ceded by another call with iopt=1 or iopt=0.
c          unchanged on exit.
c  m     : integer. on entry m must specify the number of data points.
c          m >= (kx+1)*(ky+1). unchanged on exit.
c  x     : real array of dimension at least (m).
c  y     : real array of dimension at least (m).
c  z     : real array of dimension at least (m).
c          before entry, x(i),y(i),z(i) must be set to the co-ordinates
c          of the i-th data point, for i=1,...,m. the order of the data
c          points is immaterial. unchanged on exit.
c  w     : real array of dimension at least (m). before entry, w(i) must
c          be set to the i-th value in the set of weights. the w(i) must
c          be strictly positive. unchanged on exit.
c  xb,xe : real values. on entry xb,xe,yb and ye must specify the bound-
c  yb,ye   aries of the rectangular approximation domain.
c          xb<=x(i)<=xe,yb<=y(i)<=ye,i=1,...,m. unchanged on exit.
c  kx,ky : integer values. on entry kx and ky must specify the degrees
c          of the spline. 1<=kx,ky<=5. it is recommended to use bicubic
c          (kx=ky=3) splines. unchanged on exit.
c  s     : real. on entry (in case iopt>=0) s must specify the smoothing
c          factor. s >=0. unchanged on exit.
c          for advice on the choice of s see further comments
c  nxest : integer. unchanged on exit.
c  nyest : integer. unchanged on exit.
c          on entry, nxest and nyest must specify an upper bound for the
c          number of knots required in the x- and y-directions respect.
c          these numbers will also determine the storage space needed by
c          the routine. nxest >= 2*(kx+1), nyest >= 2*(ky+1).
c          in most practical situation nxest = kx+1+sqrt(m/2), nyest =
c          ky+1+sqrt(m/2) will be sufficient. see also further comments.
c  nmax  : integer. on entry nmax must specify the actual dimension of
c          the arrays tx and ty. nmax >= nxest, nmax >=nyest.
c          unchanged on exit.
c  eps   : real.
c          on entry, eps must specify a threshold for determining the
c          effective rank of an over-determined linear system of equat-
c          ions. 0 < eps < 1.  if the number of decimal digits in the
c          computer representation of a real number is q, then 10**(-q)
c          is a suitable value for eps in most practical applications.
c          unchanged on exit.
c  nx    : integer.
c          unless ier=10 (in case iopt >=0), nx will contain the total
c          number of knots with respect to the x-variable, of the spline
c          approximation returned. if the computation mode iopt=1 is
c          used, the value of nx should be left unchanged between sub-
c          sequent calls.
c          in case iopt=-1, the value of nx should be specified on entry
c  tx    : real array of dimension nmax.
c          on succesful exit, this array will contain the knots of the
c          spline with respect to the x-variable, i.e. the position of
c          the interior knots tx(kx+2),...,tx(nx-kx-1) as well as the
c          position of the additional knots tx(1)=...=tx(kx+1)=xb and
c          tx(nx-kx)=...=tx(nx)=xe needed for the b-spline representat.
c          if the computation mode iopt=1 is used, the values of tx(1),
c          ...,tx(nx) should be left unchanged between subsequent calls.
c          if the computation mode iopt=-1 is used, the values tx(kx+2),
c          ...tx(nx-kx-1) must be supplied by the user, before entry.
c          see also the restrictions (ier=10).
c  ny    : integer.
c          unless ier=10 (in case iopt >=0), ny will contain the total
c          number of knots with respect to the y-variable, of the spline
c          approximation returned. if the computation mode iopt=1 is
c          used, the value of ny should be left unchanged between sub-
c          sequent calls.
c          in case iopt=-1, the value of ny should be specified on entry
c  ty    : real array of dimension nmax.
c          on succesful exit, this array will contain the knots of the
c          spline with respect to the y-variable, i.e. the position of
c          the interior knots ty(ky+2),...,ty(ny-ky-1) as well as the
c          position of the additional knots ty(1)=...=ty(ky+1)=yb and
c          ty(ny-ky)=...=ty(ny)=ye needed for the b-spline representat.
c          if the computation mode iopt=1 is used, the values of ty(1),
c          ...,ty(ny) should be left unchanged between subsequent calls.
c          if the computation mode iopt=-1 is used, the values ty(ky+2),
c          ...ty(ny-ky-1) must be supplied by the user, before entry.
c          see also the restrictions (ier=10).
c  c     : real array of dimension at least (nxest-kx-1)*(nyest-ky-1).
c          on succesful exit, c contains the coefficients of the spline
c          approximation s(x,y)
c  fp    : real. unless ier=10, fp contains the weighted sum of
c          squared residuals of the spline approximation returned.
c  wrk1  : real array of dimension (lwrk1). used as workspace.
c          if the computation mode iopt=1 is used the value of wrk1(1)
c          should be left unchanged between subsequent calls.
c          on exit wrk1(2),wrk1(3),...,wrk1(1+(nx-kx-1)*(ny-ky-1)) will
c          contain the values d(i)/max(d(i)),i=1,...,(nx-kx-1)*(ny-ky-1)
c          with d(i) the i-th diagonal element of the reduced triangular
c          matrix for calculating the b-spline coefficients. it includes
c          those elements whose square is less than eps,which are treat-
c          ed as 0 in the case of presumed rank deficiency (ier<-2).
c  lwrk1 : integer. on entry lwrk1 must specify the actual dimension of
c          the array wrk1 as declared in the calling (sub)program.
c          lwrk1 must not be too small. let
c            u = nxest-kx-1, v = nyest-ky-1, km = max(kx,ky)+1,
c            ne = max(nxest,nyest), bx = kx*v+ky+1, by = ky*u+kx+1,
c            if(bx.le.by) b1 = bx, b2 = b1+v-ky
c            if(bx.gt.by) b1 = by, b2 = b1+u-kx  then
c          lwrk1 >= u*v*(2+b1+b2)+2*(u+v+km*(m+ne)+ne-kx-ky)+b2+1
c  wrk2  : real array of dimension (lwrk2). used as workspace, but
c          only in the case a rank deficient system is encountered.
c  lwrk2 : integer. on entry lwrk2 must specify the actual dimension of
c          the array wrk2 as declared in the calling (sub)program.
c          lwrk2 > 0 . a save upper boundfor lwrk2 = u*v*(b2+1)+b2
c          where u,v and b2 are as above. if there are enough data
c          points, scattered uniformly over the approximation domain
c          and if the smoothing factor s is not too small, there is a
c          good chance that this extra workspace is not needed. a lot
c          of memory might therefore be saved by setting lwrk2=1.
c          (see also ier > 10)
c  iwrk  : integer array of dimension (kwrk). used as workspace.
c  kwrk  : integer. on entry kwrk must specify the actual dimension of
c          the array iwrk as declared in the calling (sub)program.
c          kwrk >= m+(nxest-2*kx-1)*(nyest-2*ky-1).
c  ier   : integer. unless the routine detects an error, ier contains a
c          non-positive value on exit, i.e.
c   ier=0  : normal return. the spline returned has a residual sum of
c            squares fp such that abs(fp-s)/s <= tol with tol a relat-
c            ive tolerance set to 0.001 by the program.
c   ier=-1 : normal return. the spline returned is an interpolating
c            spline (fp=0).
c   ier=-2 : normal return. the spline returned is the weighted least-
c            squares polynomial of degrees kx and ky. in this extreme
c            case fp gives the upper bound for the smoothing factor s.
c   ier<-2 : warning. the coefficients of the spline returned have been
c            computed as the minimal norm least-squares solution of a
c            (numerically) rank deficient system. (-ier) gives the rank.
c            especially if the rank deficiency which can be computed as
c            (nx-kx-1)*(ny-ky-1)+ier, is large the results may be inac-
c            curate. they could also seriously depend on the value of
c            eps.
c   ier=1  : error. the required storage space exceeds the available
c            storage space, as specified by the parameters nxest and
c            nyest.
c            probably causes : nxest or nyest too small. if these param-
c            eters are already large, it may also indicate that s is
c            too small
c            the approximation returned is the weighted least-squares
c            spline according to the current set of knots.
c            the parameter fp gives the corresponding weighted sum of
c            squared residuals (fp>s).
c   ier=2  : error. a theoretically impossible result was found during
c            the iteration proces for finding a smoothing spline with
c            fp = s. probably causes : s too small or badly chosen eps.
c            there is an approximation returned but the corresponding
c            weighted sum of squared residuals does not satisfy the
c            condition abs(fp-s)/s < tol.
c   ier=3  : error. the maximal number of iterations maxit (set to 20
c            by the program) allowed for finding a smoothing spline
c            with fp=s has been reached. probably causes : s too small
c            there is an approximation returned but the corresponding
c            weighted sum of squared residuals does not satisfy the
c            condition abs(fp-s)/s < tol.
c   ier=4  : error. no more knots can be added because the number of
c            b-spline coefficients (nx-kx-1)*(ny-ky-1) already exceeds
c            the number of data points m.
c            probably causes : either s or m too small.
c            the approximation returned is the weighted least-squares
c            spline according to the current set of knots.
c            the parameter fp gives the corresponding weighted sum of
c            squared residuals (fp>s).
c   ier=5  : error. no more knots can be added because the additional
c            knot would (quasi) coincide with an old one.
c            probably causes : s too small or too large a weight to an
c            inaccurate data point.
c            the approximation returned is the weighted least-squares
c            spline according to the current set of knots.
c            the parameter fp gives the corresponding weighted sum of
c            squared residuals (fp>s).
c   ier=10 : error. on entry, the input data are controlled on validity
c            the following restrictions must be satisfied.
c            -1<=iopt<=1, 1<=kx,ky<=5, m>=(kx+1)*(ky+1), nxest>=2*kx+2,
c            nyest>=2*ky+2, 0<eps<1, nmax>=nxest, nmax>=nyest,
c            xb<=x(i)<=xe, yb<=y(i)<=ye, w(i)>0, i=1,...,m
c            lwrk1 >= u*v*(2+b1+b2)+2*(u+v+km*(m+ne)+ne-kx-ky)+b2+1
c            kwrk >= m+(nxest-2*kx-1)*(nyest-2*ky-1)
c            if iopt=-1: 2*kx+2<=nx<=nxest
c                        xb<tx(kx+2)<tx(kx+3)<...<tx(nx-kx-1)<xe
c                        2*ky+2<=ny<=nyest
c                        yb<ty(ky+2)<ty(ky+3)<...<ty(ny-ky-1)<ye
c            if iopt>=0: s>=0
c            if one of these conditions is found to be violated,control
c            is immediately repassed to the calling program. in that
c            case there is no approximation returned.
c   ier>10 : error. lwrk2 is too small, i.e. there is not enough work-
c            space for computing the minimal least-squares solution of
c            a rank deficient system of linear equations. ier gives the
c            requested value for lwrk2. there is no approximation re-
c            turned but, having saved the information contained in nx,
c            ny,tx,ty,wrk1, and having adjusted the value of lwrk2 and
c            the dimension of the array wrk2 accordingly, the user can
c            continue at the point the program was left, by calling
c            surfit with iopt=1.
c
c further comments:
c  by means of the parameter s, the user can control the tradeoff
c   between closeness of fit and smoothness of fit of the approximation.
c   if s is too large, the spline will be too smooth and signal will be
c   lost ; if s is too small the spline will pick up too much noise. in
c   the extreme cases the program will return an interpolating spline if
c   s=0 and the weighted least-squares polynomial (degrees kx,ky)if s is
c   very large. between these extremes, a properly chosen s will result
c   in a good compromise between closeness of fit and smoothness of fit.
c   to decide whether an approximation, corresponding to a certain s is
c   satisfactory the user is highly recommended to inspect the fits
c   graphically.
c   recommended values for s depend on the weights w(i). if these are
c   taken as 1/d(i) with d(i) an estimate of the standard deviation of
c   z(i), a good s-value should be found in the range (m-sqrt(2*m),m+
c   sqrt(2*m)). if nothing is known about the statistical error in z(i)
c   each w(i) can be set equal to one and s determined by trial and
c   error, taking account of the comments above. the best is then to
c   start with a very large value of s ( to determine the least-squares
c   polynomial and the corresponding upper bound fp0 for s) and then to
c   progressively decrease the value of s ( say by a factor 10 in the
c   beginning, i.e. s=fp0/10, fp0/100,...and more carefully as the
c   approximation shows more detail) to obtain closer fits.
c   to choose s very small is strongly discouraged. this considerably
c   increases computation time and memory requirements. it may also
c   cause rank-deficiency (ier<-2) and endager numerical stability.
c   to economize the search for a good s-value the program provides with
c   different modes of computation. at the first call of the routine, or
c   whenever he wants to restart with the initial set of knots the user
c   must set iopt=0.
c   if iopt=1 the program will continue with the set of knots found at
c   the last call of the routine. this will save a lot of computation
c   time if surfit is called repeatedly for different values of s.
c   the number of knots of the spline returned and their location will
c   depend on the value of s and on the complexity of the shape of the
c   function underlying the data. if the computation mode iopt=1
c   is used, the knots returned may also depend on the s-values at
c   previous calls (if these were smaller). therefore, if after a number
c   of trials with different s-values and iopt=1, the user can finally
c   accept a fit as satisfactory, it may be worthwhile for him to call
c   surfit once more with the selected value for s but now with iopt=0.
c   indeed, surfit may then return an approximation of the same quality
c   of fit but with fewer knots and therefore better if data reduction
c   is also an important objective for the user.
c   the number of knots may also depend on the upper bounds nxest and
c   nyest. indeed, if at a certain stage in surfit the number of knots
c   in one direction (say nx) has reached the value of its upper bound
c   (nxest), then from that moment on all subsequent knots are added
c   in the other (y) direction. this may indicate that the value of
c   nxest is too small. on the other hand, it gives the user the option
c   of limiting the number of knots the routine locates in any direction
c   for example, by setting nxest=2*kx+2 (the lowest allowable value for
c   nxest), the user can indicate that he wants an approximation which
c   is a simple polynomial of degree kx in the variable x.
c
c  other subroutines required:
c    fpback,fpbspl,fpsurf,fpdisc,fpgivs,fprank,fprati,fprota,fporde
c
c  references:
c   dierckx p. : an algorithm for surface fitting with spline functions
c                ima j. numer. anal. 1 (1981) 267-283.
c   dierckx p. : an algorithm for surface fitting with spline functions
c                report tw50, dept. computer science,k.u.leuven, 1980.
c   dierckx p. : curve and surface fitting with splines, monographs on
c                numerical analysis, oxford university press, 1993.
c
c  author:
c    p.dierckx
c    dept. computer science, k.u. leuven
c    celestijnenlaan 200a, b-3001 heverlee, belgium.
c    e-mail : Paul.Dierckx@@cs.kuleuven.ac.be
c
c  creation date : may 1979
c  latest update : march 1987
c
c  ..
c  ..scalar arguments..
      real*8 xb,xe,yb,ye,s,eps,fp
      integer iopt,m,kx,ky,nxest,nyest,nmax,nx,ny,lwrk1,lwrk2,kwrk,ier
c  ..array arguments..
      real*8 x(m),y(m),z(m),w(m),tx(nmax),ty(nmax),
     * c((nxest-kx-1)*(nyest-ky-1)),wrk1(lwrk1),wrk2(lwrk2)
      integer iwrk(kwrk)
c  ..local scalars..
      real*8 tol
      integer i,ib1,ib3,jb1,ki,kmax,km1,km2,kn,kwest,kx1,ky1,la,lbx,
     * lby,lco,lf,lff,lfp,lh,lq,lsx,lsy,lwest,maxit,ncest,nest,nek,
     * nminx,nminy,nmx,nmy,nreg,nrint,nxk,nyk
c  ..function references..
      integer max0
c  ..subroutine references..
c    fpsurf
c  ..
c  we set up the parameters tol and maxit.
      maxit = 20
      tol = 0.1e-02
c  before starting computations a data check is made. if the input data
c  are invalid,control is immediately repassed to the calling program.
      ier = 10
      if(eps.le.0. .or. eps.ge.1.) go to 70
      if(kx.le.0 .or. kx.gt.5) go to 70
      kx1 = kx+1
      if(ky.le.0 .or. ky.gt.5) go to 70
      ky1 = ky+1
      kmax = max0(kx,ky)
      km1 = kmax+1
      km2 = km1+1
      if(iopt.lt.(-1) .or. iopt.gt.1) go to 70
      if(m.lt.(kx1*ky1)) go to 70
      nminx = 2*kx1
      if(nxest.lt.nminx .or. nxest.gt.nmax) go to 70
      nminy = 2*ky1
      if(nyest.lt.nminy .or. nyest.gt.nmax) go to 70
      nest = max0(nxest,nyest)
      nxk = nxest-kx1
      nyk = nyest-ky1
      ncest = nxk*nyk
      nmx = nxest-nminx+1
      nmy = nyest-nminy+1
      nrint = nmx+nmy
      nreg = nmx*nmy
      ib1 = kx*nyk+ky1
      jb1 = ky*nxk+kx1
      ib3 = kx1*nyk+1
      if(ib1.le.jb1) go to 10
      ib1 = jb1
      ib3 = ky1*nxk+1
  10  lwest = ncest*(2+ib1+ib3)+2*(nrint+nest*km2+m*km1)+ib3
      kwest = m+nreg
      if(lwrk1.lt.lwest .or. kwrk.lt.kwest) go to 70
      if(xb.ge.xe .or. yb.ge.ye) go to 70
      do 20 i=1,m
        if(w(i).le.0.) go to 70
        if(x(i).lt.xb .or. x(i).gt.xe) go to 70
        if(y(i).lt.yb .or. y(i).gt.ye) go to 70
  20  continue
      if(iopt.ge.0) go to 50
      if(nx.lt.nminx .or. nx.gt.nxest) go to 70
      nxk = nx-kx1
      tx(kx1) = xb
      tx(nxk+1) = xe
      do 30 i=kx1,nxk
        if(tx(i+1).le.tx(i)) go to 70
  30  continue
      if(ny.lt.nminy .or. ny.gt.nyest) go to 70
      nyk = ny-ky1
      ty(ky1) = yb
      ty(nyk+1) = ye
      do 40 i=ky1,nyk
        if(ty(i+1).le.ty(i)) go to 70
  40  continue
      go to 60
  50  if(s.lt.0.) go to 70
  60  ier = 0
c  we partition the working space and determine the spline approximation
      kn = 1
      ki = kn+m
      lq = 2
      la = lq+ncest*ib3
      lf = la+ncest*ib1
      lff = lf+ncest
      lfp = lff+ncest
      lco = lfp+nrint
      lh = lco+nrint
      lbx = lh+ib3
      nek = nest*km2
      lby = lbx+nek
      lsx = lby+nek
      lsy = lsx+m*km1
      call fpsurf(iopt,m,x,y,z,w,xb,xe,yb,ye,kx,ky,s,nxest,nyest,
     * eps,tol,maxit,nest,km1,km2,ib1,ib3,ncest,nrint,nreg,nx,tx,
     * ny,ty,c,fp,wrk1(1),wrk1(lfp),wrk1(lco),wrk1(lf),wrk1(lff),
     * wrk1(la),wrk1(lq),wrk1(lbx),wrk1(lby),wrk1(lsx),wrk1(lsy),
     * wrk1(lh),iwrk(ki),iwrk(kn),wrk2,lwrk2,ier)
  70  return
      end
      subroutine fpback(a,z,n,k,c,nest)
c  subroutine fpback calculates the solution of the system of
c  equations a*c = z with a a n x n upper triangular matrix
c  of bandwidth k.
c  ..
c  ..scalar arguments..
      integer n,k,nest
c  ..array arguments..
      real*8 a(nest,k),z(n),c(n)
c  ..local scalars..
      real*8 store
      integer i,i1,j,k1,l,m
c  ..
      k1 = k-1
      c(n) = z(n)/a(n,1)
      i = n-1
      if(i.eq.0) go to 30
      do 20 j=2,n
        store = z(i)
        i1 = k1
        if(j.le.k1) i1 = j-1
        m = i
        do 10 l=1,i1
          m = m+1
          store = store-c(m)*a(i,l+1)
  10    continue
        c(i) = store/a(i,1)
        i = i-1
  20  continue
  30  return
      end
      subroutine fpsurf(iopt,m,x,y,z,w,xb,xe,yb,ye,kxx,kyy,s,nxest,
     * nyest,eta,tol,maxit,nmax,km1,km2,ib1,ib3,nc,intest,nrest,
     * nx0,tx,ny0,ty,c,fp,fp0,fpint,coord,f,ff,a,q,bx,by,spx,spy,h,
     * index,nummer,wrk,lwrk,ier)
c  ..
c  ..scalar arguments..
      real*8 xb,xe,yb,ye,s,eta,tol,fp,fp0
      integer iopt,m,kxx,kyy,nxest,nyest,maxit,nmax,km1,km2,ib1,ib3,
     * nc,intest,nrest,nx0,ny0,lwrk,ier
c  ..array arguments..
      real*8 x(m),y(m),z(m),w(m),tx(nmax),ty(nmax),c(nc),fpint(intest),
     * coord(intest),f(nc),ff(nc),a(nc,ib1),q(nc,ib3),bx(nmax,km2),
     * by(nmax,km2),spx(m,km1),spy(m,km1),h(ib3),wrk(lwrk)
      integer index(nrest),nummer(m)
c  ..local scalars..
      real*8 acc,arg,cos,dmax,fac1,fac2,fpmax,fpms,f1,f2,f3,hxi,p,pinv,
     * piv,p1,p2,p3,sigma,sin,sq,store,wi,x0,x1,y0,y1,zi,eps,
     * rn,one,con1,con9,con4,half,ten
      integer i,iband,iband1,iband3,iband4,ibb,ichang,ich1,ich3,ii,
     * in,irot,iter,i1,i2,i3,j,jrot,jxy,j1,kx,kx1,kx2,ky,ky1,ky2,l,
     * la,lf,lh,lwest,lx,ly,l1,l2,n,ncof,nk1x,nk1y,nminx,nminy,nreg,
     * nrint,num,num1,nx,nxe,nxx,ny,nye,nyy,n1,rank
c  ..local arrays..
      real*8 hx(6),hy(6)
c  ..function references..
      real*8 abs,fprati,sqrt
      integer min0
c  ..subroutine references..
c    fpback,fpbspl,fpgivs,fpdisc,fporde,fprank,fprota
c  ..
c  set constants
      one = 0.1e+01
      con1 = 0.1e0
      con9 = 0.9e0
      con4 = 0.4e-01
      half = 0.5e0
      ten = 0.1e+02
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c part 1: determination of the number of knots and their position.     c
c ****************************************************************     c
c given a set of knots we compute the least-squares spline sinf(x,y),  c
c and the corresponding weighted sum of squared residuals fp=f(p=inf). c
c if iopt=-1  sinf(x,y) is the requested approximation.                c
c if iopt=0 or iopt=1 we check whether we can accept the knots:        c
c   if fp <=s we will continue with the current set of knots.          c
c   if fp > s we will increase the number of knots and compute the     c
c      corresponding least-squares spline until finally  fp<=s.        c
c the initial choice of knots depends on the value of s and iopt.      c
c   if iopt=0 we first compute the least-squares polynomial of degree  c
c     kx in x and ky in y; nx=nminx=2*kx+2 and ny=nminy=2*ky+2.        c
c     fp0=f(0) denotes the corresponding weighted sum of squared       c
c     residuals                                                        c
c   if iopt=1 we start with the knots found at the last call of the    c
c     routine, except for the case that s>=fp0; then we can compute    c
c     the least-squares polynomial directly.                           c
c eventually the independent variables x and y (and the corresponding  c
c parameters) will be switched if this can reduce the bandwidth of the c
c system to be solved.                                                 c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  ichang denotes whether(1) or not(-1) the directions have been inter-
c  changed.
      ichang = -1
      x0 = xb
      x1 = xe
      y0 = yb
      y1 = ye
      kx = kxx
      ky = kyy
      kx1 = kx+1
      ky1 = ky+1
      nxe = nxest
      nye = nyest
      eps = sqrt(eta)
      if(iopt.lt.0) go to 20
c  calculation of acc, the absolute tolerance for the root of f(p)=s.
      acc = tol*s
      if(iopt.eq.0) go to 10
      if(fp0.gt.s) go to 20
c  initialization for the least-squares polynomial.
  10  nminx = 2*kx1
      nminy = 2*ky1
      nx = nminx
      ny = nminy
      ier = -2
      go to 30
  20  nx = nx0
      ny = ny0
c  main loop for the different sets of knots. m is a save upper bound
c  for the number of trials.
  30  do 420 iter=1,m
c  find the position of the additional knots which are needed for the
c  b-spline representation of s(x,y).
        l = nx
        do 40 i=1,kx1
          tx(i) = x0
          tx(l) = x1
          l = l-1
  40    continue
        l = ny
        do 50 i=1,ky1
          ty(i) = y0
          ty(l) = y1
          l = l-1
  50    continue
c  find nrint, the total number of knot intervals and nreg, the number
c  of panels in which the approximation domain is subdivided by the
c  intersection of knots.
        nxx = nx-2*kx1+1
        nyy = ny-2*ky1+1
        nrint = nxx+nyy
        nreg = nxx*nyy
c  find the bandwidth of the observation matrix a.
c  if necessary, interchange the variables x and y, in order to obtain
c  a minimal bandwidth.
        iband1 = kx*(ny-ky1)+ky
        l = ky*(nx-kx1)+kx
        if(iband1.le.l) go to 130
        iband1 = l
        ichang = -ichang
        do 60 i=1,m
          store = x(i)
          x(i) = y(i)
          y(i) = store
  60    continue
        store = x0
        x0 = y0
        y0 = store
        store = x1
        x1 = y1
        y1 = store
        n = min0(nx,ny)
        do 70 i=1,n
          store = tx(i)
          tx(i) = ty(i)
          ty(i) = store
  70    continue
        n1 = n+1
        if(nx-ny) 80,120,100
  80    do 90 i=n1,ny
          tx(i) = ty(i)
  90    continue
        go to 120
 100    do 110 i=n1,nx
          ty(i) = tx(i)
 110    continue
 120    l = nx
        nx = ny
        ny = l
        l = nxe
        nxe = nye
        nye = l
        l = nxx
        nxx = nyy
        nyy = l
        l = kx
        kx = ky
        ky = l
        kx1 = kx+1
        ky1 = ky+1
 130    iband = iband1+1
c  arrange the data points according to the panel they belong to.
        call fporde(x,y,m,kx,ky,tx,nx,ty,ny,nummer,index,nreg)
c  find ncof, the number of b-spline coefficients.
        nk1x = nx-kx1
        nk1y = ny-ky1
        ncof = nk1x*nk1y
c  initialize the observation matrix a.
        do 140 i=1,ncof
          f(i) = 0.
          do 140 j=1,iband
            a(i,j) = 0.
 140    continue
c  initialize the sum of squared residuals.
        fp = 0.
c  fetch the data points in the new order. main loop for the
c  different panels.
        do 250 num=1,nreg
c  fix certain constants for the current panel; jrot records the column
c  number of the first non-zero element in a row of the observation
c  matrix according to a data point of the panel.
          num1 = num-1
          lx = num1/nyy
          l1 = lx+kx1
          ly = num1-lx*nyy
          l2 = ly+ky1
          jrot = lx*nk1y+ly
c  test whether there are still data points in the panel.
          in = index(num)
 150      if(in.eq.0) go to 250
c  fetch a new data point.
          wi = w(in)
          zi = z(in)*wi
c  evaluate for the x-direction, the (kx+1) non-zero b-splines at x(in).
          call fpbspl(tx,nx,kx,x(in),l1,hx)
c  evaluate for the y-direction, the (ky+1) non-zero b-splines at y(in).
          call fpbspl(ty,ny,ky,y(in),l2,hy)
c  store the value of these b-splines in spx and spy respectively.
          do 160 i=1,kx1
            spx(in,i) = hx(i)
 160      continue
          do 170 i=1,ky1
            spy(in,i) = hy(i)
 170      continue
c  initialize the new row of observation matrix.
          do 180 i=1,iband
            h(i) = 0.
 180      continue
c  calculate the non-zero elements of the new row by making the cross
c  products of the non-zero b-splines in x- and y-direction.
          i1 = 0
          do 200 i=1,kx1
            hxi = hx(i)
            j1 = i1
            do 190 j=1,ky1
              j1 = j1+1
              h(j1) = hxi*hy(j)*wi
 190        continue
            i1 = i1+nk1y
 200      continue
c  rotate the row into triangle by givens transformations .
          irot = jrot
          do 220 i=1,iband
            irot = irot+1
            piv = h(i)
            if(piv.eq.0.) go to 220
c  calculate the parameters of the givens transformation.
            call fpgivs(piv,a(irot,1),cos,sin)
c  apply that transformation to the right hand side.
            call fprota(cos,sin,zi,f(irot))
            if(i.eq.iband) go to 230
c  apply that transformation to the left hand side.
            i2 = 1
            i3 = i+1
            do 210 j=i3,iband
              i2 = i2+1
              call fprota(cos,sin,h(j),a(irot,i2))
 210        continue
 220      continue
c  add the contribution of the row to the sum of squares of residual
c  right hand sides.
 230      fp = fp+zi**2
c  find the number of the next data point in the panel.
 240      in = nummer(in)
          go to 150
 250    continue
c  find dmax, the maximum value for the diagonal elements in the reduced
c  triangle.
        dmax = 0.
        do 260 i=1,ncof
          if(a(i,1).le.dmax) go to 260
          dmax = a(i,1)
 260    continue
c  check whether the observation matrix is rank deficient.
        sigma = eps*dmax
        do 270 i=1,ncof
          if(a(i,1).le.sigma) go to 280
 270    continue
c  backward substitution in case of full rank.
        call fpback(a,f,ncof,iband,c,nc)
        rank = ncof
        do 275 i=1,ncof
          q(i,1) = a(i,1)/dmax
 275    continue
        go to 300
c  in case of rank deficiency, find the minimum norm solution.
c  check whether there is sufficient working space
 280    lwest = ncof*iband+ncof+iband
        if(lwrk.lt.lwest) go to 780
        do 290 i=1,ncof
          ff(i) = f(i)
          do 290 j=1,iband
            q(i,j) = a(i,j)
 290    continue
        lf =1
        lh = lf+ncof
        la = lh+iband
        call fprank(q,ff,ncof,iband,nc,sigma,c,sq,rank,wrk(la),
     *    wrk(lf),wrk(lh))
        do 295 i=1,ncof
          q(i,1) = q(i,1)/dmax
 295    continue
c  add to the sum of squared residuals, the contribution of reducing
c  the rank.
        fp = fp+sq
 300    if(ier.eq.(-2)) fp0 = fp
c  test whether the least-squares spline is an acceptable solution.
        if(iopt.lt.0) go to 820
        fpms = fp-s
        if(abs(fpms).le.acc) if(fp) 815,815,820
c  test whether we can accept the choice of knots.
        if(fpms.lt.0.) go to 430
c  test whether we cannot further increase the number of knots.
        if(ncof.gt.m) go to 790
        ier = 0
c  search where to add a new knot.
c  find for each interval the sum of squared residuals fpint for the
c  data points having the coordinate belonging to that knot interval.
c  calculate also coord which is the same sum, weighted by the position
c  of the data points considered.
 310    do 320 i=1,nrint
          fpint(i) = 0.
          coord(i) = 0.
 320    continue
        do 360 num=1,nreg
          num1 = num-1
          lx = num1/nyy
          l1 = lx+1
          ly = num1-lx*nyy
          l2 = ly+1+nxx
          jrot = lx*nk1y+ly
          in = index(num)
 330      if(in.eq.0) go to 360
          store = 0.
          i1 = jrot
          do 350 i=1,kx1
            hxi = spx(in,i)
            j1 = i1
            do 340 j=1,ky1
              j1 = j1+1
              store = store+hxi*spy(in,j)*c(j1)
 340        continue
            i1 = i1+nk1y
 350      continue
          store = (w(in)*(z(in)-store))**2
          fpint(l1) = fpint(l1)+store
          coord(l1) = coord(l1)+store*x(in)
          fpint(l2) = fpint(l2)+store
          coord(l2) = coord(l2)+store*y(in)
          in = nummer(in)
          go to 330
 360    continue
c  find the interval for which fpint is maximal on the condition that
c  there still can be added a knot.
 370    l = 0
        fpmax = 0.
        l1 = 1
        l2 = nrint
        if(nx.eq.nxe) l1 = nxx+1
        if(ny.eq.nye) l2 = nxx
        if(l1.gt.l2) go to 810
        do 380 i=l1,l2
          if(fpmax.ge.fpint(i)) go to 380
          l = i
          fpmax = fpint(i)
 380    continue
c  test whether we cannot further increase the number of knots.
        if(l.eq.0) go to 785
c  calculate the position of the new knot.
        arg = coord(l)/fpint(l)
c  test in what direction the new knot is going to be added.
        if(l.gt.nxx) go to 400
c  addition in the x-direction.
        jxy = l+kx1
        fpint(l) = 0.
        fac1 = tx(jxy)-arg
        fac2 = arg-tx(jxy-1)
        if(fac1.gt.(ten*fac2) .or. fac2.gt.(ten*fac1)) go to 370
        j = nx
        do 390 i=jxy,nx
          tx(j+1) = tx(j)
          j = j-1
 390    continue
        tx(jxy) = arg
        nx = nx+1
        go to 420
c  addition in the y-direction.
 400    jxy = l+ky1-nxx
        fpint(l) = 0.
        fac1 = ty(jxy)-arg
        fac2 = arg-ty(jxy-1)
        if(fac1.gt.(ten*fac2) .or. fac2.gt.(ten*fac1)) go to 370
        j = ny
        do 410 i=jxy,ny
          ty(j+1) = ty(j)
          j = j-1
 410    continue
        ty(jxy) = arg
        ny = ny+1
c  restart the computations with the new set of knots.
 420  continue
c  test whether the least-squares polynomial is a solution of our
c  approximation problem.
 430  if(ier.eq.(-2)) go to 830
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c part 2: determination of the smoothing spline sp(x,y)                c
c *****************************************************                c
c we have determined the number of knots and their position. we now    c
c compute the b-spline coefficients of the smoothing spline sp(x,y).   c
c the observation matrix a is extended by the rows of a matrix,        c
c expressing that sp(x,y) must be a polynomial of degree kx in x and   c
c ky in y. the corresponding weights of these additional rows are set  c
c to 1./p.  iteratively we than have to determine the value of p       c
c such that f(p)=sum((w(i)*(z(i)-sp(x(i),y(i))))**2) be = s.           c
c we already know that the least-squares polynomial corresponds to     c
c p=0  and that the least-squares spline corresponds to p=infinity.    c
c the iteration process which is proposed here makes use of rational   c
c interpolation. since f(p) is a convex and strictly decreasing        c
c function of p, it can be approximated by a rational function r(p)=   c
c (u*p+v)/(p+w). three values of p(p1,p2,p3) with corresponding values c
c of f(p) (f1=f(p1)-s,f2=f(p2)-s,f3=f(p3)-s) are used to calculate the c
c new value of p such that r(p)=s. convergence is guaranteed by taking c
c f1 > 0 and f3 < 0.                                                   c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      kx2 = kx1+1
c  test whether there are interior knots in the x-direction.
      if(nk1x.eq.kx1) go to 440
c  evaluate the discotinuity jumps of the kx-th order derivative of
c  the b-splines at the knots tx(l),l=kx+2,...,nx-kx-1.
      call fpdisc(tx,nx,kx2,bx,nmax)
 440  ky2 = ky1 + 1
c  test whether there are interior knots in the y-direction.
      if(nk1y.eq.ky1) go to 450
c  evaluate the discontinuity jumps of the ky-th order derivative of
c  the b-splines at the knots ty(l),l=ky+2,...,ny-ky-1.
      call fpdisc(ty,ny,ky2,by,nmax)
c  initial value for p.
 450  p1 = 0.
      f1 = fp0-s
      p3 = -one
      f3 = fpms
      p = 0.
      do 460 i=1,ncof
        p = p+a(i,1)
 460  continue
      rn = ncof
      p = rn/p
c  find the bandwidth of the extended observation matrix.
      iband3 = kx1*nk1y
      iband4 = iband3 +1
      ich1 = 0
      ich3 = 0
c  iteration process to find the root of f(p)=s.
      do 770 iter=1,maxit
        pinv = one/p
c  store the triangularized observation matrix into q.
        do 480 i=1,ncof
          ff(i) = f(i)
          do 470 j=1,iband
            q(i,j) = a(i,j)
 470      continue
          ibb = iband+1
          do 480 j=ibb,iband4
            q(i,j) = 0.
 480    continue
        if(nk1y.eq.ky1) go to 560
c  extend the observation matrix with the rows of a matrix, expressing
c  that for x=cst. sp(x,y) must be a polynomial in y of degree ky.
        do 550 i=ky2,nk1y
          ii = i-ky1
          do 550 j=1,nk1x
c  initialize the new row.
            do 490 l=1,iband
              h(l) = 0.
 490        continue
c  fill in the non-zero elements of the row. jrot records the column
c  number of the first non-zero element in the row.
            do 500 l=1,ky2
              h(l) = by(ii,l)*pinv
 500        continue
            zi = 0.
            jrot = (j-1)*nk1y+ii
c  rotate the new row into triangle by givens transformations without
c  square roots.
            do 540 irot=jrot,ncof
              piv = h(1)
              i2 = min0(iband1,ncof-irot)
              if(piv.eq.0.) if(i2) 550,550,520
c  calculate the parameters of the givens transformation.
              call fpgivs(piv,q(irot,1),cos,sin)
c  apply that givens transformation to the right hand side.
              call fprota(cos,sin,zi,ff(irot))
              if(i2.eq.0) go to 550
c  apply that givens transformation to the left hand side.
              do 510 l=1,i2
                l1 = l+1
                call fprota(cos,sin,h(l1),q(irot,l1))
 510          continue
 520          do 530 l=1,i2
                h(l) = h(l+1)
 530          continue
              h(i2+1) = 0.
 540        continue
 550    continue
 560    if(nk1x.eq.kx1) go to 640
c  extend the observation matrix with the rows of a matrix expressing
c  that for y=cst. sp(x,y) must be a polynomial in x of degree kx.
        do 630 i=kx2,nk1x
          ii = i-kx1
          do 630 j=1,nk1y
c  initialize the new row
            do 570 l=1,iband4
              h(l) = 0.
 570        continue
c  fill in the non-zero elements of the row. jrot records the column
c  number of the first non-zero element in the row.
            j1 = 1
            do 580 l=1,kx2
              h(j1) = bx(ii,l)*pinv
              j1 = j1+nk1y
 580        continue
            zi = 0.
            jrot = (i-kx2)*nk1y+j
c  rotate the new row into triangle by givens transformations .
            do 620 irot=jrot,ncof
              piv = h(1)
              i2 = min0(iband3,ncof-irot)
              if(piv.eq.0.) if(i2) 630,630,600
c  calculate the parameters of the givens transformation.
              call fpgivs(piv,q(irot,1),cos,sin)
c  apply that givens transformation to the right hand side.
              call fprota(cos,sin,zi,ff(irot))
              if(i2.eq.0) go to 630
c  apply that givens transformation to the left hand side.
              do 590 l=1,i2
                l1 = l+1
                call fprota(cos,sin,h(l1),q(irot,l1))
 590          continue
 600          do 610 l=1,i2
                h(l) = h(l+1)
 610          continue
              h(i2+1) = 0.
 620        continue
 630    continue
c  find dmax, the maximum value for the diagonal elements in the
c  reduced triangle.
 640    dmax = 0.
        do 650 i=1,ncof
          if(q(i,1).le.dmax) go to 650
          dmax = q(i,1)
 650    continue
c  check whether the matrix is rank deficient.
        sigma = eps*dmax
        do 660 i=1,ncof
          if(q(i,1).le.sigma) go to 670
 660    continue
c  backward substitution in case of full rank.
        call fpback(q,ff,ncof,iband4,c,nc)
        rank = ncof
        go to 675
c  in case of rank deficiency, find the minimum norm solution.
 670    lwest = ncof*iband4+ncof+iband4
        if(lwrk.lt.lwest) go to 780
        lf = 1
        lh = lf+ncof
        la = lh+iband4
        call fprank(q,ff,ncof,iband4,nc,sigma,c,sq,rank,wrk(la),
     *   wrk(lf),wrk(lh))
 675    do 680 i=1,ncof
          q(i,1) = q(i,1)/dmax
 680    continue
c  compute f(p).
        fp = 0.
        do 720 num = 1,nreg
          num1 = num-1
          lx = num1/nyy
          ly = num1-lx*nyy
          jrot = lx*nk1y+ly
          in = index(num)
 690      if(in.eq.0) go to 720
          store = 0.
          i1 = jrot
          do 710 i=1,kx1
            hxi = spx(in,i)
            j1 = i1
            do 700 j=1,ky1
              j1 = j1+1
              store = store+hxi*spy(in,j)*c(j1)
 700        continue
            i1 = i1+nk1y
 710      continue
          fp = fp+(w(in)*(z(in)-store))**2
          in = nummer(in)
          go to 690
 720    continue
c  test whether the approximation sp(x,y) is an acceptable solution.
        fpms = fp-s
        if(abs(fpms).le.acc) go to 820
c  test whether the maximum allowable number of iterations has been
c  reached.
        if(iter.eq.maxit) go to 795
c  carry out one more step of the iteration process.
        p2 = p
        f2 = fpms
        if(ich3.ne.0) go to 740
        if((f2-f3).gt.acc) go to 730
c  our initial choice of p is too large.
        p3 = p2
        f3 = f2
        p = p*con4
        if(p.le.p1) p = p1*con9 + p2*con1
        go to 770
 730    if(f2.lt.0.) ich3 = 1
 740    if(ich1.ne.0) go to 760
        if((f1-f2).gt.acc) go to 750
c  our initial choice of p is too small
        p1 = p2
        f1 = f2
        p = p/con4
        if(p3.lt.0.) go to 770
        if(p.ge.p3) p = p2*con1 + p3*con9
        go to 770
 750    if(f2.gt.0.) ich1 = 1
c  test whether the iteration process proceeds as theoretically
c  expected.
 760    if(f2.ge.f1 .or. f2.le.f3) go to 800
c  find the new value of p.
        p = fprati(p1,f1,p2,f2,p3,f3)
 770  continue
c  error codes and messages.
 780  ier = lwest
      go to 830
 785  ier = 5
      go to 830
 790  ier = 4
      go to 830
 795  ier = 3
      go to 830
 800  ier = 2
      go to 830
 810  ier = 1
      go to 830
 815  ier = -1
      fp = 0.
 820  if(ncof.ne.rank) ier = -rank
c  test whether x and y are in the original order.
 830  if(ichang.lt.0) go to 930
c  if not, interchange x and y once more.
      l1 = 1
      do 840 i=1,nk1x
        l2 = i
        do 840 j=1,nk1y
          f(l2) = c(l1)
          l1 = l1+1
          l2 = l2+nk1x
 840  continue
      do 850 i=1,ncof
        c(i) = f(i)
 850  continue
      do 860 i=1,m
        store = x(i)
        x(i) = y(i)
        y(i) = store
 860  continue
      n = min0(nx,ny)
      do 870 i=1,n
        store = tx(i)
        tx(i) = ty(i)
        ty(i) = store
 870  continue
      n1 = n+1
      if(nx-ny) 880,920,900
 880  do 890 i=n1,ny
        tx(i) = ty(i)
 890  continue
      go to 920
 900  do 910 i=n1,nx
        ty(i) = tx(i)
 910  continue
 920  l = nx
      nx = ny
      ny = l
 930  if(iopt.lt.0) go to 940
      nx0 = nx
      ny0 = ny
 940  return
      end

      subroutine fpdisc(t,n,k2,b,nest)
c  subroutine fpdisc calculates the discontinuity jumps of the kth
c  derivative of the b-splines of degree k at the knots t(k+2)..t(n-k-1)
c  ..scalar arguments..
      integer n,k2,nest
c  ..array arguments..
      real*8 t(n),b(nest,k2)
c  ..local scalars..
      real*8 an,fac,prod
      integer i,ik,j,jk,k,k1,l,lj,lk,lmk,lp,nk1,nrint
c  ..local array..
      real*8 h(12)
c  ..
      k1 = k2-1
      k = k1-1
      nk1 = n-k1
      nrint = nk1-k
      an = nrint
      fac = an/(t(nk1+1)-t(k1))
      do 40 l=k2,nk1
        lmk = l-k1
        do 10 j=1,k1
          ik = j+k1
          lj = l+j
          lk = lj-k2
          h(j) = t(l)-t(lk)
          h(ik) = t(l)-t(lj)
  10    continue
        lp = lmk
        do 30 j=1,k2
          jk = j
          prod = h(j)
          do 20 i=1,k
            jk = jk+1
            prod = prod*h(jk)*fac
  20      continue
          lk = lp+k1
          b(lmk,j) = (t(lk)-t(lp))/prod
          lp = lp+1
  30    continue
  40  continue
      return
      end
      subroutine fpgivs(piv,ww,cos,sin)
c  subroutine fpgivs calculates the parameters of a givens
c  transformation .
c  ..
c  ..scalar arguments..
      real*8 piv,ww,cos,sin
c  ..local scalars..
      real*8 dd,one,store
c  ..function references..
      real*8 abs,sqrt
c  ..
      one = 0.1e+01
      store = abs(piv)
      if(store.ge.ww) dd = store*sqrt(one+(ww/piv)**2)
      if(store.lt.ww) dd = ww*sqrt(one+(piv/ww)**2)
      cos = ww/dd
      sin = piv/dd
      ww = dd
      return
      end
      subroutine fprank(a,f,n,m,na,tol,c,sq,rank,aa,ff,h)
c  subroutine fprank finds the minimum norm solution of a least-
c  squares problem in case of rank deficiency.
c
c  input parameters:
c    a : array, which contains the non-zero elements of the observation
c        matrix after triangularization by givens transformations.
c    f : array, which contains the transformed right hand side.
c    n : integer,wich contains the dimension of a.
c    m : integer, which denotes the bandwidth of a.
c  tol : real value, giving a threshold to determine the rank of a.
c
c  output parameters:
c    c : array, which contains the minimum norm solution.
c   sq : real value, giving the contribution of reducing the rank
c        to the sum of squared residuals.
c rank : integer, which contains the rank of matrix a.
c
c  ..scalar arguments..
      integer n,m,na,rank
      real*8 tol,sq
c  ..array arguments..
      real*8 a(na,m),f(n),c(n),aa(n,m),ff(n),h(m)
c  ..local scalars..
      integer i,ii,ij,i1,i2,j,jj,j1,j2,j3,k,kk,m1,nl
      real*8 cos,fac,piv,sin,yi
      double precision store,stor1,stor2,stor3
c  ..function references..
      integer min0
c  ..subroutine references..
c    fpgivs,fprota
c  ..
      m1 = m-1
c  the rank deficiency nl is considered to be the number of sufficient
c  small diagonal elements of a.
      nl = 0
      sq = 0.
      do 90 i=1,n
        if(a(i,1).gt.tol) go to 90
c  if a sufficient small diagonal element is found, we put it to
c  zero. the remainder of the row corresponding to that zero diagonal
c  element is then rotated into triangle by givens rotations .
c  the rank deficiency is increased by one.
        nl = nl+1
        if(i.eq.n) go to 90
        yi = f(i)
        do 10 j=1,m1
          h(j) = a(i,j+1)
  10    continue
        h(m) = 0.
        i1 = i+1
        do 60 ii=i1,n
          i2 = min0(n-ii,m1)
          piv = h(1)
          if(piv.eq.0.) go to 30
          call fpgivs(piv,a(ii,1),cos,sin)
          call fprota(cos,sin,yi,f(ii))
          if(i2.eq.0) go to 70
          do 20 j=1,i2
            j1 = j+1
            call fprota(cos,sin,h(j1),a(ii,j1))
            h(j) = h(j1)
  20      continue
          go to 50
  30      if(i2.eq.0) go to 70
          do 40 j=1,i2
            h(j) = h(j+1)
  40      continue
  50      h(i2+1) = 0.
  60    continue
c  add to the sum of squared residuals the contribution of deleting
c  the row with small diagonal element.
  70    sq = sq+yi**2
  90  continue
c  rank denotes the rank of a.
      rank = n-nl
c  let b denote the (rank*n) upper trapezoidal matrix which can be
c  obtained from the (n*n) upper triangular matrix a by deleting
c  the rows and interchanging the columns corresponding to a zero
c  diagonal element. if this matrix is factorized using givens
c  transformations as  b = (r) (u)  where
c    r is a (rank*rank) upper triangular matrix,
c    u is a (rank*n) orthonormal matrix
c  then the minimal least-squares solution c is given by c = b' v,
c  where v is the solution of the system  (r) (r)' v = g  and
c  g denotes the vector obtained from the old right hand side f, by
c  removing the elements corresponding to a zero diagonal element of a.
c  initialization.
      do 100 i=1,rank
        do 100 j=1,m
          aa(i,j) = 0.
 100  continue
c  form in aa the upper triangular matrix obtained from a by
c  removing rows and columns with zero diagonal elements. form in ff
c  the new right hand side by removing the elements of the old right
c  hand side corresponding to a deleted row.
      ii = 0
      do 120 i=1,n
        if(a(i,1).le.tol) go to 120
        ii = ii+1
        ff(ii) = f(i)
        aa(ii,1) = a(i,1)
        jj = ii
        kk = 1
        j = i
        j1 = min0(j-1,m1)
        if(j1.eq.0) go to 120
        do 110 k=1,j1
          j = j-1
          if(a(j,1).le.tol) go to 110
          kk = kk+1
          jj = jj-1
          aa(jj,kk) = a(j,k+1)
 110    continue
 120  continue
c  form successively in h the columns of a with a zero diagonal element.
      ii = 0
      do 200 i=1,n
        ii = ii+1
        if(a(i,1).gt.tol) go to 200
        ii = ii-1
        if(ii.eq.0) go to 200
        jj = 1
        j = i
        j1 = min0(j-1,m1)
        do 130 k=1,j1
          j = j-1
          if(a(j,1).le.tol) go to 130
          h(jj) = a(j,k+1)
          jj = jj+1
 130    continue
        do 140 kk=jj,m
          h(kk) = 0.
 140    continue
c  rotate this column into aa by givens transformations.
        jj = ii
        do 190 i1=1,ii
          j1 = min0(jj-1,m1)
          piv = h(1)
          if(piv.ne.0.) go to 160
          if(j1.eq.0) go to 200
          do 150 j2=1,j1
            j3 = j2+1
            h(j2) = h(j3)
 150      continue
          go to 180
 160      call fpgivs(piv,aa(jj,1),cos,sin)
          if(j1.eq.0) go to 200
          kk = jj
          do 170 j2=1,j1
            j3 = j2+1
            kk = kk-1
            call fprota(cos,sin,h(j3),aa(kk,j3))
            h(j2) = h(j3)
 170      continue
 180      jj = jj-1
          h(j3) = 0.
 190    continue
 200  continue
c  solve the system (aa) (f1) = ff
      ff(rank) = ff(rank)/aa(rank,1)
      i = rank-1
      if(i.eq.0) go to 230
      do 220 j=2,rank
        store = ff(i)
        i1 = min0(j-1,m1)
        k = i
        do 210 ii=1,i1
          k = k+1
          stor1 = ff(k)
          stor2 = aa(i,ii+1)
          store = store-stor1*stor2
 210    continue
        stor1 = aa(i,1)
        ff(i) = store/stor1
        i = i-1
 220  continue
c  solve the system  (aa)' (f2) = f1
 230  ff(1) = ff(1)/aa(1,1)
      if(rank.eq.1) go to 260
      do 250 j=2,rank
        store = ff(j)
        i1 = min0(j-1,m1)
        k = j
        do 240 ii=1,i1
          k = k-1
          stor1 = ff(k)
          stor2 = aa(k,ii+1)
          store = store-stor1*stor2
 240    continue
        stor1 = aa(j,1)
        ff(j) = store/stor1
 250  continue
c  premultiply f2 by the transpoze of a.
 260  k = 0
      do 280 i=1,n
        store = 0.
        if(a(i,1).gt.tol) k = k+1
        j1 = min0(i,m)
        kk = k
        ij = i+1
        do 270 j=1,j1
          ij = ij-1
          if(a(ij,1).le.tol) go to 270
          stor1 = a(ij,j)
          stor2 = ff(kk)
          store = store+stor1*stor2
          kk = kk-1
 270    continue
        c(i) = store
 280  continue
c  add to the sum of squared residuals the contribution of putting
c  to zero the small diagonal elements of matrix (a).
      stor3 = 0.
      do 310 i=1,n
        if(a(i,1).gt.tol) go to 310
        store = f(i)
        i1 = min0(n-i,m1)
        if(i1.eq.0) go to 300
        do 290 j=1,i1
          ij = i+j
          stor1 = c(ij)
          stor2 = a(i,j+1)
          store = store-stor1*stor2
 290    continue
 300    fac = a(i,1)*c(i)
        stor1 = a(i,1)
        stor2 = c(i)
        stor1 = stor1*stor2
        stor3 = stor3+stor1*(stor1-store-store)
 310  continue
      fac = stor3
      sq = sq+fac
      return
      end

      real*8 function fprati(p1,f1,p2,f2,p3,f3)
c  given three points (p1,f1),(p2,f2) and (p3,f3), function fprati
c  gives the value of p such that the rational interpolating function
c  of the form r(p) = (u*p+v)/(p+w) equals zero at p.
c  ..
c  ..scalar arguments..
      real*8 p1,f1,p2,f2,p3,f3
c  ..local scalars..
      real*8 h1,h2,h3,p
c  ..
      if(p3.gt.0.) go to 10
c  value of p in case p3 = infinity.
      p = (p1*(f1-f3)*f2-p2*(f2-f3)*f1)/((f1-f2)*f3)
      go to 20
c  value of p in case p3 ^= infinity.
  10  h1 = f1*(f2-f3)
      h2 = f2*(f3-f1)
      h3 = f3*(f1-f2)
      p = -(p1*p2*h3+p2*p3*h1+p3*p1*h2)/(p1*h1+p2*h2+p3*h3)
c  adjust the value of p1,f1,p3 and f3 such that f1 > 0 and f3 < 0.
  20  if(f2.lt.0.) go to 30
      p1 = p2
      f1 = f2
      go to 40
  30  p3 = p2
      f3 = f2
  40  fprati = p
      return
      end
      subroutine fprota(cos,sin,a,b)
c  subroutine fprota applies a givens rotation to a and b.
c  ..
c  ..scalar arguments..
      real*8 cos,sin,a,b
c ..local scalars..
      real*8 stor1,stor2
c  ..
      stor1 = a
      stor2 = b
      b = cos*stor2+sin*stor1
      a = cos*stor1-sin*stor2
      return
      end
      subroutine fporde(x,y,m,kx,ky,tx,nx,ty,ny,nummer,index,nreg)
c  subroutine fporde sorts the data points (x(i),y(i)),i=1,2,...,m
c  according to the panel tx(l)<=x<tx(l+1),ty(k)<=y<ty(k+1), they belong
c  to. for each panel a stack is constructed  containing the numbers
c  of data points lying inside; index(j),j=1,2,...,nreg points to the
c  first data point in the jth panel while nummer(i),i=1,2,...,m gives
c  the number of the next data point in the panel.
c  ..
c  ..scalar arguments..
      integer m,kx,ky,nx,ny,nreg
c  ..array arguments..
      real*8 x(m),y(m),tx(nx),ty(ny)
      integer nummer(m),index(nreg)
c  ..local scalars..
      real*8 xi,yi
      integer i,im,k,kx1,ky1,k1,l,l1,nk1x,nk1y,num,nyy
c  ..
      kx1 = kx+1
      ky1 = ky+1
      nk1x = nx-kx1
      nk1y = ny-ky1
      nyy = nk1y-ky
      do 10 i=1,nreg
        index(i) = 0
  10  continue
      do 60 im=1,m
        xi = x(im)
        yi = y(im)
        l = kx1
        l1 = l+1
  20    if(xi.lt.tx(l1) .or. l.eq.nk1x) go to 30
        l = l1
        l1 = l+1
        go to 20
  30    k = ky1
        k1 = k+1
  40    if(yi.lt.ty(k1) .or. k.eq.nk1y) go to 50
        k = k1
        k1 = k+1
        go to 40
  50    num = (l-kx1)*nyy+k-ky
        nummer(im) = index(num)
        index(num) = im
  60  continue
      return
      end
@}

The c-wrapper routines:

@o splines.c
@{
/* These were Copyrighted by Enthought Inc */
/* They  seemed  to  give permission  for  me to copy out the 
   parts I wanted. I promise that if it goes wrong it was not
   their fault.  They are nice people and not responsible for
   errors I might have introduced                  - JPW 2004  */

#include <Python.h>
#include "Numeric/arrayobject.h"  

#define BISPEV bispev_
#define PARDER parder_
#define SURFIT surfit_

void BISPEV(double*,int*,double*,int*,double*,int*,int*,double*,int*,
                      double*,int*,double*,double*,int*,int*,int*,int*);
void PARDER(double*,int*,double*,int*,double*,int*,int*,int*,int*,double*,
                  int*,double*,int*,double*,double*,int*,int*,int*,int*);
void SURFIT(int*,int*,double*,double*,double*,double*,double*,double*,double*,double*,
            int*,int*,double*,int*,int*,int*,double*,int*,double*,int*,double*,double*,
            double*,double*,int*,double*,int*,int*,int*,int*);

static char doc_bispev[] = " [z,ier] = _bispev(tx,ty,c,kx,ky,x,y,nux,nuy)";

static PyObject *_bispev(PyObject *dummy, PyObject *args) {
  int nx,ny,kx,ky,mx,my,lwrk,*iwrk,kwrk,ier,lwa,mxy,nux,nuy;
  double *tx,*ty,*c,*x,*y,*z,*wrk,*wa = NULL;
  PyArrayObject *ap_x = NULL,*ap_y = NULL,*ap_z = NULL,*ap_tx = NULL,\
    *ap_ty = NULL,*ap_c = NULL;
  PyObject *x_py = NULL,*y_py = NULL,*c_py = NULL,*tx_py = NULL,*ty_py = NULL;
  if (!PyArg_ParseTuple(args, "OOOiiOOii",&tx_py,&ty_py,&c_py,&kx,&ky,
                        &x_py,&y_py,&nux,&nuy))
    return NULL;
  ap_x = (PyArrayObject *)PyArray_ContiguousFromObject(x_py, PyArray_DOUBLE, 0, 1);
  ap_y = (PyArrayObject *)PyArray_ContiguousFromObject(y_py, PyArray_DOUBLE, 0, 1);
  ap_c = (PyArrayObject *)PyArray_ContiguousFromObject(c_py, PyArray_DOUBLE, 0, 1);
  ap_tx = (PyArrayObject *)PyArray_ContiguousFromObject(tx_py, PyArray_DOUBLE, 0, 1);
  ap_ty = (PyArrayObject *)PyArray_ContiguousFromObject(ty_py, PyArray_DOUBLE, 0, 1);
  if (ap_x == NULL || ap_y == NULL || ap_c == NULL || ap_tx == NULL \
      || ap_ty == NULL) goto fail;  
  x = (double *) ap_x->data;
  y = (double *) ap_y->data;
  c = (double *) ap_c->data;
  tx = (double *) ap_tx->data;
  ty = (double *) ap_ty->data;
  nx = ap_tx->dimensions[0];
  ny = ap_ty->dimensions[0];
  mx = ap_x->dimensions[0];
  my = ap_y->dimensions[0];
  mxy = mx*my;
  ap_z = (PyArrayObject *)PyArray_FromDims(1,&mxy,PyArray_DOUBLE);
  z = (double *) ap_z->data;
  if (nux || nuy) 
    lwrk = mx*(kx+1-nux)+my*(ky+1-nuy)+(nx-kx-1)*(ny-ky-1);
  else
    lwrk = mx*(kx+1)+my*(ky+1);
  kwrk = mx+my;
  lwa = lwrk+kwrk;
  if ((wa = (double *)malloc(lwa*sizeof(double)))==NULL) {
    PyErr_NoMemory();
    goto fail;
  }
  wrk = wa;
  iwrk = (int *)(wrk+lwrk);
  if (nux || nuy)
    PARDER(tx,&nx,ty,&ny,c,&kx,&ky,&nux,&nuy,x,&mx,y,&my,z,wrk,&lwrk,iwrk,&kwrk,&ier);
  else
    BISPEV(tx,&nx,ty,&ny,c,&kx,&ky,x,&mx,y,&my,z,wrk,&lwrk,iwrk,&kwrk,&ier);
  if (wa) free(wa);
  Py_DECREF(ap_x);
  Py_DECREF(ap_y);
  Py_DECREF(ap_c);
  Py_DECREF(ap_tx);
  Py_DECREF(ap_ty);
  return Py_BuildValue("Ni",PyArray_Return(ap_z),ier);
  fail:
  if (wa) free(wa);
  Py_XDECREF(ap_x);
  Py_XDECREF(ap_y);
  Py_XDECREF(ap_z);
  Py_XDECREF(ap_c);
  Py_XDECREF(ap_tx);
  Py_XDECREF(ap_ty);
  return NULL;
}

static char doc_surfit[] =
 " [tx,ty,c,o] = "
 "_surfit(x,y,z,w,xb,xe,yb,ye,kx,ky,iopt,s,eps,tx,ty,nxest,nyest,wrk,lwrk1,lwrk2)";

static PyObject *_surfit(PyObject *dummy, PyObject *args) {
  int iopt,m,kx,ky,nxest,nyest,nx,ny,lwrk1,lwrk2,*iwrk,kwrk,ier,lwa,nxo,nyo,\
    i,lc,lcest,nmax;
  double *x,*y,*z,*w,xb,xe,yb,ye,s,*tx,*ty,*c,fp,*wrk1,*wrk2,*wa = NULL,eps;
  PyArrayObject *ap_x = NULL,*ap_y = NULL,*ap_z,*ap_w = NULL,\
    *ap_tx = NULL,*ap_ty = NULL,*ap_c = NULL;
  PyArrayObject *ap_wrk = NULL;
  PyObject *x_py = NULL,*y_py = NULL,*z_py = NULL,*w_py = NULL,\
    *tx_py = NULL,*ty_py = NULL;
  PyObject *wrk_py=NULL;
  nx=ny=ier=nxo=nyo=0;
  if (!PyArg_ParseTuple(args, "OOOOddddiiiddOOiiOii",\
			&x_py,&y_py,&z_py,&w_py,&xb,&xe,\
			&yb,&ye,&kx,&ky,&iopt,&s,&eps,&tx_py,&ty_py,&nxest,&nyest,\
			&wrk_py,&lwrk1,&lwrk2)) return NULL;
  ap_x = (PyArrayObject *)PyArray_ContiguousFromObject(x_py, PyArray_DOUBLE, 0, 1);
  ap_y = (PyArrayObject *)PyArray_ContiguousFromObject(y_py, PyArray_DOUBLE, 0, 1);
  ap_z = (PyArrayObject *)PyArray_ContiguousFromObject(z_py, PyArray_DOUBLE, 0, 1);
  ap_w = (PyArrayObject *)PyArray_ContiguousFromObject(w_py, PyArray_DOUBLE, 0, 1);
  ap_wrk=(PyArrayObject *)PyArray_ContiguousFromObject(wrk_py, PyArray_DOUBLE, 0, 1);
  /*ap_iwrk=(PyArrayObject *)PyArray_ContiguousFromObject(iwrk_py, PyArray_INT, 0, 1);*/
  if (ap_x == NULL || ap_y == NULL || ap_z == NULL || ap_w == NULL \
      || ap_wrk == NULL) goto fail;
  x = (double *) ap_x->data;
  y = (double *) ap_y->data;
  z = (double *) ap_z->data;
  w = (double *) ap_w->data;
  m = ap_x->dimensions[0];
  nmax=nxest;
  if (nmax<nyest) nmax=nyest;
  lcest=(nxest-kx-1)*(nyest-ky-1);
  kwrk=m+(nxest-2*kx-1)*(nyest-2*ky-1);
  lwa = 2*nmax+lcest+lwrk1+lwrk2+kwrk;
  if ((wa = (double *)malloc(lwa*sizeof(double)))==NULL) {
    PyErr_NoMemory();
    goto fail;
  }
  tx = wa;
  ty = tx + nmax;
  c = ty + nmax;
  wrk1 = c + lcest;
  iwrk = (int *)(wrk1 + lwrk1);
  wrk2 = (double *)(iwrk+kwrk);
  if (iopt) {
    ap_tx=(PyArrayObject *)PyArray_ContiguousFromObject(tx_py, PyArray_DOUBLE, 0, 1);
    ap_ty=(PyArrayObject *)PyArray_ContiguousFromObject(ty_py, PyArray_DOUBLE, 0, 1);
    if (ap_tx == NULL || ap_ty == NULL) goto fail;
    nx = nxo = ap_tx->dimensions[0];
    ny = nyo = ap_ty->dimensions[0];
    memcpy(tx,ap_tx->data,nx*sizeof(double));
    memcpy(ty,ap_ty->data,ny*sizeof(double));
  }
  if (iopt==1) {
    lc = (nx-kx-1)*(ny-ky-1);
    memcpy(wrk1,ap_wrk->data,lc*sizeof(double));
    /*memcpy(iwrk,ap_iwrk->data,n*sizeof(int));*/
  }
  SURFIT(&iopt,&m,x,y,z,w,&xb,&xe,&yb,&ye,&kx,&ky,&s,&nxest,&nyest,&nmax,
        &eps,&nx,tx,&ny,ty,c,&fp,wrk1,&lwrk1,wrk2,&lwrk2,iwrk,&kwrk,&ier);
  i=0;
  while ((ier>10) && (i++<5)) {
    lwrk2=ier;
    if ((wrk2 = (double *)malloc(lwrk2*sizeof(double)))==NULL) {
      PyErr_NoMemory();
      goto fail;
    }
    SURFIT(&iopt,&m,x,y,z,w,&xb,&xe,&yb,&ye,&kx,&ky,&s,&nxest,&nyest,&nmax,
          &eps,&nx,tx,&ny,ty,c,&fp,wrk1,&lwrk1,wrk2,&lwrk2,iwrk,&kwrk,&ier);
    if (wrk2) free(wrk2);
  }
  if (ier==10) goto fail;
  lc = (nx-kx-1)*(ny-ky-1);
  ap_tx = (PyArrayObject *)PyArray_FromDims(1,&nx,PyArray_DOUBLE);
  ap_ty = (PyArrayObject *)PyArray_FromDims(1,&ny,PyArray_DOUBLE);
  ap_c = (PyArrayObject *)PyArray_FromDims(1,&lc,PyArray_DOUBLE);
  if (ap_tx == NULL || ap_ty == NULL || ap_c == NULL) goto fail;
  if ((iopt==0)||(nx>nxo)||(ny>nyo)) {
    ap_wrk = (PyArrayObject *)PyArray_FromDims(1,&lc,PyArray_DOUBLE);
    if (ap_wrk == NULL) goto fail;
    /*ap_iwrk = (PyArrayObject *)PyArray_FromDims(1,&n,PyArray_INT);*/
  }
  memcpy(ap_tx->data,tx,nx*sizeof(double));
  memcpy(ap_ty->data,ty,ny*sizeof(double));
  memcpy(ap_c->data,c,lc*sizeof(double));
  memcpy(ap_wrk->data,wrk1,lc*sizeof(double));
  /*memcpy(ap_iwrk->data,iwrk,n*sizeof(int));*/
  if (wa) free(wa);
  Py_DECREF(ap_x);
  Py_DECREF(ap_y);
  Py_DECREF(ap_z);
  Py_DECREF(ap_w);
  return Py_BuildValue("NNN{s:N,s:i,s:d}",PyArray_Return(ap_tx),\
		       PyArray_Return(ap_ty),PyArray_Return(ap_c),\
		       "wrk",PyArray_Return(ap_wrk),\
		       "ier",ier,"fp",fp);
  fail:
  if (wa) free(wa);
  Py_XDECREF(ap_x);
  Py_XDECREF(ap_y);
  Py_XDECREF(ap_z);
  Py_XDECREF(ap_w);
  Py_XDECREF(ap_tx);
  Py_XDECREF(ap_ty);
  Py_XDECREF(ap_wrk);
  /*Py_XDECREF(ap_iwrk);*/
  return NULL;
}



static struct PyMethodDef splines_methods[] = {
	{"_bispev", _bispev, METH_VARARGS, doc_bispev},
   {"_surfit", _surfit, METH_VARARGS, doc_surfit},

	{NULL,          NULL, 0, NULL}
	};

void init_splines(void) {
  PyObject *m, *d, *s;
  m = Py_InitModule("_splines", splines_methods);
  import_array();
  d = PyModule_GetDict(m);
  s = PyString_FromString(" 0.1 ");
  PyDict_SetItemString(d, "__version__", s);
  Py_DECREF(s);
  if (PyErr_Occurred())
    Py_FatalError("can't initialize module splines");
}
@}




\subsubsection{Generating spline files}

The routine for fitting a spline to a set of points is also available.
To use this we would take a list of fitted grid positions, ask for an ideal
peak (one to start the delta\_x and delta\_y as zero) and ask for the pixel
and mask dimensions and grid rotation angle. Then when each peak is assigned
a true hole position (based on the mask geometry and pixel size) we can calculate
the offset between the measured grid position and calculated position if the 
pixels were perfectly layed out with no distortion. This set of points define
a deformation surface which we want to get the spline for. The true pixel positions 
are then the x,y from i,j array indices plus the spline function.

Need to work through an example, from a list of grid positions. Things to do are
figure out the mask spacing (in pixels) and the rotation angle of the grid. 
In principle these can be automated to give a closest fitting set of positions
before fitting the spline to get the final correction. Vignetting will be a pain,
if it is done at all. Ought to find a pincushion function, if appropriate, for the 
image intensifier.

Assuming we have a set of peak positions, perhaps from imagepro, we can try to go through
and work out the spline file stuff....

\begin{verbatim}
#Total Count:	4117

 Obj#	     Area	   Aspect	Area (polygon)	Center-X (mass)	Center-Y (mass)

    1	       34	 1.745368	 20.16667	 17.58042	 1.836231
    2	       43	 1.601456	 29.79167	 37.62301	 3.978285
    5	       32	 1.284832	    20.50	 151.6240	 2.379340
    7	       21	 1.843058	 8.430550	 275.0284	 1.048654
...
\end{verbatim}

First we read in an Imagepro "cnt" file:

@d readipcnt
@{
def readipcnt(name):
   """
   Reads an imagepro list of objects found into an array of:
      ar[0,:]=Center-X
      ar[1,:]=Center-Y
   """
   f=open(name)
   l=f.readline() # we just don't care how many blobs there were
   l=f.readline() # due to the joys of lists
   l=f.readline()
   titles=l.split("\t")
   for i in range(len(titles)):
      if titles[i].find("Center-X") != -1: xcol=i
      if titles[i].find("Center-Y") != -1: ycol=i
   l=f.readline()
   data=f.readlines()
   f.close()
   x=[]
   y=[]
   for l in data:
      x.append(float(l.split()[xcol]))
      y.append(float(l.split()[ycol]))
   return array((x,y),Float32)
@}

\subsection{Generating an imagepro style blob file}

Assuming we have run the threshold blob search (see a later section)
we should have a data image and
an image of integer peak assignments. This could be done much faster within
the c-code which does the blob searching, but for now we will make something
quick and dirty in python.

@O blobproperties.py
@{
import _connectedpixels

from Numeric import *

@< readipcnt @>

def connectedpixels(ar,t,Centreofmass=0): 
   if Centreofmass==0: return _connectedpixels._connectedpixels(ar,t)
   else : return _connectedpixels._connectedpixels(ar,t,Centreofmass=Centreofmass)

def blobcounts(data,threshold):
   blobarray= _connectedpixels._connectedpixels(data,threshold,Centreofmass=6)
   # npt, sumI, sumX, sumY, sumXI, sumYI
   # so we have size, intensity, sumX/npt=avg x, 
   blobarray[:,2]=blobarray[:,2]/blobarray[:,0]
   blobarray[:,3]=blobarray[:,3]/blobarray[:,0]
   blobarray[:,4]=blobarray[:,4]/blobarray[:,1]
   blobarray[:,5]=blobarray[:,5]/blobarray[:,1]
   return blobarray


def blobproperties(data,index):
   nblobs=maximum.reduce(ravel(index))
   area=zeros(nblobs+1,Int)
   cx=zeros(nblobs+1,Float)
   cy=zeros(nblobs+1,Float)
   for i in range(index.shape[0]):
      print "Processing row ",i,"\r",
      for j in range(index.shape[1]):
         blob=index[i,j]
         area[blob]=area[blob]+1
         cx[blob]=cx[blob]+i
         cy[blob]=cy[blob]+j
   print
   cx=cx/area
   cy=cy/area
   return area,cx,cy



if __name__=="__main__":
   import sys
   if __name__=="__main__":
      try:
         filename=sys.argv[1]
         threshold=float(sys.argv[2])
      except:
         print "Usage: %s filename threshold"%(sys.argv[0])
         sys.exit()
   if len(sys.argv)>3: outfile=open(sys.argv[3],"w")
   else: outfile=sys.stdout
   from time import time
   import testedffile
   import testbruker
   t1=time()
   if filename[-3:]=='edf': 
      head,image=testedffile.readedffile(filename)
   else:
      head,image=testbruker.readbruker(filename)
   t2=time()
   print "Time reading file =", t2-t1
   bp=blobcounts(image,threshold)
   t3=time()
   print "Blob counting took",t3-t2
   outfile.write("#Total Count: %d\n"%(bp.shape[0]-1)) 
           # skip the zeroth order pixel
   outfile.write(
       "\nObj#\tIntensity\t#pixels\tCenter-X (mass)\tCenter-Y (mass)\n\n")
   for i in range(1,bp.shape[0]):
      outfile.write("%6d\t%12.2f\t%d\t%7f\t%.7f\n"%
                    (i,int(bp[i,1]),int(bp[i,0]),bp[i,4],bp[i,5]))
   if outfile!=sys.stdout: outfile.close()
@}

\subsection{Indexing the peaks from a grid image}

We can read an array of \code{xy[2,npoints]} where there are x and y values
for each of the peaks from the cnt file. 
We now need to assign the each peak a set of x and y co-ordinates
such that each peak also has a corresponding pair of integers attached to it.

Proof that you can write fortran in any language..., eventually this should
be put into a c extension, but it is something which runs infrequently
and usually goes wrong, so that won't happen till it works really really well
and someone has more free time.

@O indexingpeaks.py
@{from Numeric import *

big=99999999          # be afraid, be very afraid


def smallest(ar,n):
   """
   Returns the positions of the n smallest elements in array
   Should be a Numeric thing to do this, it is only faster than
   sorting for very few numbers from very big arrays
   """
   ret=zeros(n,Int)
   m=minimum(ar.shape[0],n)
   if(ar.shape[0] > 3*n*n):   # if less than 3n^2
      x=maximum.reduce(ar)    # max in input
      for i in range(m):      
         ret[i]=argmin(ar)    # get min and modify ar so we don't get it next time
         ar[ret[i]]=x+ar[ret[i]]
      for i in range(m):      # put the right numbers back
         ar[ret[i]]=ar[ret[i]]-x
   else:
      t=argsort(ar)           # more elegant, but slower
      ret[:m]=t[:m]   
   return ret

def findneighbours(xy,searchrange=-1):
   """
   Find the neighbour peaks and compute distances to them
   Returns an array of nearest neighbours (9 including itself),
   sorted by distance in n[0:9,i] where i are the neighbours
   for a particular peak. Also returns an array of distances
   which correspond to the peaks in n.
   """
   r=searchrange         # search range variable
   n=zeros( (9,xy.shape[1]) ,Int)      # result neighbours
   dist=zeros( (9,xy.shape[1]), Float) # result distances
   for p in range(xy.shape[1]):
      x=xy[0,p] # x - coord
      y=xy[1,p] # y - coord
      # Compute distance to each, sort, and take the lowest values
      # was a bit slow, now using adaptive search range
      # and smallest routine to try to go faster than sort
      if r < 1:
         d=(xy[0,:]-x)**2+(xy[1,:]-y)**2   # compute distances to all peaks
         n[:,p]=smallest(d,9)  # for large arrays and few values this is faster
         dist[1:,p]=sqrt(take(d,n[1:,p]))
      else: # Hard to believe we can be so inefficient, use a search range...   
         i=array(range(xy.shape[1]))
         u=compress( xy[1,:]         < y+r , i )  # u now hold indices below y+r
         u=compress( take(xy[1,:],u) > y-r , u )  # u now has indices above y-r
         u=compress( take(xy[0,:],u) < x+r , u )  # now only where x < x+r
         u=compress( take(xy[0,:],u) > x-r , u )  # now only where x > x-r
         t=take(xy,u,1)                           # x,y data values
         d=(t[0,:]-x)**2+(t[1,:]-y)**2            # distance (squared)
         a=smallest(d,9)                          # pick out the small distances
         n[:,p]=take(u,a[0:9])                    # what if d has fewer than 9 peaks?
         dist[1:,p]=sqrt(take(d,a[1:9]))
      if (p+1)%100 == 0: # p is the peak being treated
         ad=add.reduce(ravel(dist[1:,p-99:p]))/(8.*99)
         s="Found neighbours for %d peaks, average distance=%f pixels\r" % (p,ad)
         print s,
         # Use average distance *3 for adaptive search range
         r=ad*3  
   print "\nTreated ",p,"peaks OK!"
   return n,dist

def assignneighbours(p,n,xy,ind,ad,todo):
   """
   Given an assigned peak (p), a neighbour array (n) and the imagepro
   array xy this tries to fill in "ind", the array of indices.
   ad is the average peak to peak distance and todo is a list
   of peaks whose neighbours must be examined
   """
   if ind[0,p]!=big:        # big implies peak is unindexed
      for i in range(1,5):
         q=n[i,p] # Which peak
         if ind[0,q]!=big: continue # peak is already assigned
         dx=(xy[0,p]-xy[0,q])/ad
         dy=(xy[1,p]-xy[1,q])/ad
         if fabs(dy)<0.2 and dx > 0.4 and dx < 1.2:
            ind[0,q]=ind[0,p]+1
            ind[1,q]=ind[1,p]
         elif fabs(dy)<0.2 and dx < -0.4 and dx > -1.2:
            ind[0,q]=ind[0,p]-1
            ind[1,q]=ind[1,p]
         elif fabs(dx)<0.2 and dy > 0.4 and dy < 1.2:
            ind[0,q]=ind[0,p]
            ind[1,q]=ind[1,p]+1
         elif fabs(dx)<0.2 and dy < -0.4 and dy > -1.2:
            ind[0,q]=ind[0,p]
            ind[1,q]=ind[1,p]-1
         if ind[0,q]!=big: 
            todo.append(q)

def getims(name):
   """
   Reads an imagepro distortion/grid/object search file and tries to assign the peaks
   ready for a spline.
   """
   xy=readipcnt(name)          # Imagepro's points
   n,d=findneighbours(xy)
   ind=zeros(xy.shape,Int)+big # Will hold the i,j peak indexing 
   p=xy.shape[1]/2             # middle peak
   ind[0,p]=0 
   ind[1,p]=0
   ad=add.reduce(ravel(d[1:,:]))/(8.*d.shape[1])
   todo=[]
   todo.append(p)
   i=0
   while len(todo)>0:
      next=todo.pop(0)
      assignneighbours(next,n,xy,ind,ad,todo)
      if (i+1)%100 == 0:
         s="Assigning peak number %d\r" % i
         print s,
      i=i+1
   print "\n"
   print sum(where(ind[0,:]==big,1,0)),"peaks not assigned"
   xm=minimum.reduce(where(ind[0,:]!=big,ind[0,:],big))
   ym=minimum.reduce(where(ind[0,:]!=big,ind[1,:],big))
   xr=maximum.reduce(where(ind[0,:]!=big,ind[0,:],0))-xm+1
   yr=maximum.reduce(where(ind[0,:]!=big,ind[1,:],0))-ym+1
   xpeakarray=zeros((xr,yr),Float32)
   ypeakarray=zeros((xr,yr),Float32)
   for i in range(ind.shape[1]):
      if ind[0,i]!=big:
         xpeakarray[ind[0,i]-xm,ind[1,i]-ym]=xy[0,i]
         ypeakarray[ind[0,i]-xm,ind[1,i]-ym]=xy[1,i]
   return xpeakarray,ypeakarray,xy,ind,n

# TO DO  - alter and fix peak indexing if it goes wrong
#        - let the user delete unwanted peaks

@< readipcnt @>

if __name__=="__main__":
   import sys

   xpeakarray,ypeakarray,xy,ind,n=getims(sys.argv[1])
   print xpeakarray.shape
   print ypeakarray.shape
   
@}

Only the \code{findneighbours} part of this routine takes a long
time, and I guess that might end up being invarient enough to
warrant doing in C. You'd just have a local list of the nine smallest 
found so far, and the biggest of the nine smallest peaks for testing against.
Think about it. 

\subsection{Examining and modifying the results of the indexing}

We have a pair of arrays, called xpeak array and ypeak array which
we would like to draw on the screen. They respectively contain the 
x and y positions of the spots which have been indexed. We also have
an index array which contains the indices of each spot, and a flag
value called big which indicates unindexed spots.

Drawing lines on a canvas which connect the indexed and adjacent points
would be a useful way to show what has happened. We need to know
the i,j indexing of each point and it's co-ordinates. If we use the
neighbour array this task becomes much easier - we just check
each of the neighbours for each peak to see if they are 
indexed and then draw the lines accordingly.


\subsection{Fitting a spline file} 


@D calcideal
@{
def calcideal(theta,siz,xp,yp):
   """
   calculates an ideal set of grid positions for an array
   with rotation angle starting at theta with grid spacing starting at siz
   xp and yp are the arrays of observed peak positions
   """
   from LinearAlgebra import inverse
   keepgoing=1; i=0
   f=lambda i,j: i*st+j*ct
   g=lambda i,j: i*ct+j*st
   # compute derivate of positions w.r.t angle
   df=lambda i,j:  i*ct-j*st
   dg=lambda i,j: -i*st+j*ct
   np=multiply.reduce(xp.shape)
   # obs - calc
   b=zeros((2*np),Float)     # obs-calc
   d=zeros((2,2*np),Float)   # derivatives
   while keepgoing==1:
      ct=cos(theta)
      st=sin(theta) 
      sh=xp.shape
      xpi=array(fromfunction(f,sh),Float)*siz
      ypi=array(fromfunction(g,sh),Float)*siz
      dx=add.reduce(ravel(where(xp>0,xp-xpi,0)))
      dy=add.reduce(ravel(where(yp>0,yp-ypi,0)))
      npoint=add.reduce(ravel(where(xp>0,1,0)))
      # Centre the calculated grid on top of the observed one:
      xpi=xpi+float(dx/npoint)
      ypi=ypi+float(dy/npoint)   
      # difference
      b[:np]=ravel(where(xp>0,xp-xpi,0))
      b[np:]=ravel(where(yp>0,yp-ypi,0))
      # compute derivative of difference w.r.t size
      d[0,:np]=ravel(where(xpi>0,-xpi/siz,0))
      d[0,np:]=ravel(where(ypi>0,-ypi/siz,0))
      d[1,:np]=-ravel(array(fromfunction(df,sh),Float)*siz)
      d[1,np:]=-ravel(array(fromfunction(dg,sh),Float)*siz)
      # zero out derivatives for unobserved points
      d[1,:np]=where(ravel(xp)>0,d[1,:np],0)
      d[1,np:]=where(ravel(yp)>0,d[1,np:],0)
      x=matrixmultiply(d,b) # 2*1
      a=matrixmultiply(d,transpose(d)) # 2*2
      x=matrixmultiply(inverse(a),x)
      siz=siz+x[0]
      theta=theta+x[1]
      print i, siz, inverse(a)[0,0], theta, inverse(a)[1,1]
      if (x[1] < inverse(a)[1,1]) and i>3:
         print "Converged"
         keepgoing=0
      if i>10:
         keepgoing=0
      i=i+1
   return xpi,ypi,d,b

xp,yp,xy,ind=getims("gridpoints.cnt")
theta=0.01
siz=16
xpi,ypi,d,b=calcideal(0.001,15,xp,yp)
# Optimising the grid spacing and rotation angle.....
# the obs positions are returned from getims and we can calculate
# some from calcideal, if we guess at initial theta and size
#
# where(obs>0,(obs-calc),0.)   gives the differences
# where obs is both the x and y distortion images
# Least squares needs....
# A.x=b   
# A=d(calc)/d(parameter)
# x=to be found
# b=differences


@}







\chapter{Dealing with raw data - peak searching and fitting}

In handling the initial raw data coming for various experiments
it is not unreasonable to expect a computer program to be
able to have a quick look and to find the blindingly obvious
peaks at least. These will then be available for getting
initial values for unit cell parameters, diffractometer 
calibrations and peak shape functions. The question 
"what numbers should I use to start my refinement?" ought not
to arise in simple cases at least.

This chapter should eventually provide a program for 
integrating powder data into correlated integrated intensities. 

\section{1D peak searching}

There is something in J.Appl.Cryst to look up about using 
a Savitsky-Golay filter and examining the second derivatives?
FIXME

We can also estimate the background and ask for a threshold level
above the background.

Given an arrays of data in x and y we go through the y array with some
window size estimating the background as the smallest value within
that window. This value is then subtracted off of the y data to 
give a new array. Then some threshold (ideally some number of sigma)
is used to determine if a datapoint is high enough to be a peak. 
Connected data points are then assumed to be a part of the same 
peak and the centre and width can be estimated from the points
above the threshold level.

\section{2D peak searching}

The algorithm of the previous section can be extended to a 2 dimensional
image by using a background determined from a 2D area and a threshold 
again supplied by a user.
Before going onto the actual peak search algorithm we can take a brief
detour on the world of disjoint sets, which are an extraordinarily 
useful thing for assigning pixels to peaks quickly.

\subsection{Disjoint sets}

Sadly I could not find the library function for making a disjoint set
data structure. From "Introduction to Algorithms" (Cormen, Leiserson,
Rivesa and Stein):

\begin{verbatim}
MAKE-SET(x)
1 p[x] <- x
2 rank[x] <- 0

UNION(x,y)
1 LINK(FIND-SET(x),FIND-SET(y))

LINK(x,y)
1 if(rank[x]>rank[y])
2   then p[y] <- x
3   else p[x] <- y
4      if rank[x]==rank[y]
5         then rank[y] <- rank[y]+1

FIND-SET(x)
1 if x != p[x]
2    then p[x] <- FIND-SET(p[x])
3 return p[x]
\end{verbatim}

First a python implementation - make a list of x, p[x] values using a dictionary.
I'm sure this can be done with arrays in a better way as we always
index the dictionary as p[x=integer]. Also the recursive FIND-SET function
can be replaced with a loop for languages which do not allow recursion 
(eg fortran). Something which traverses down the tree until it finds 
p[x]==x and then returns back up setting all the p[x] equal to the 
top of the tree. Or alternatively travels down it a second time once
the answer is know.

eg"a=0; x0=0; while p[x] != x:{x=p[x];} a=x; x=x0; while p[x]!=x: { x=p[x]; p[x]=a}"

@D disjointset
@{
class disjointset:
   """Implement a disjoint set data structure for 
   fast pixel to peak assignments
   """
   def __init__(self):
      self.p={}     # initialise empty dictionaries
#     self.rank={}
      self.F={}
      self.length=0
   def makeset(self):
      self.length=self.length+1
      x=self.length # always monotonically increasing
      self.p[x]=x
#     self.rank[x]=0
      return x
   def makeunion(self,x,y):
      x=int(x)
      y=int(y)
      self.link(self.findset(x),self.findset(y))
   def link(self,x,y):
     x=int(x)
     y=int(y)
     if x>y: self.p[x]=y
     if x<y: self.p[y]=x
#     if(self.rank[x]>self.rank[y]):
#        self.p[y] = x
#     else:
#        self.p[x] = y
#        if self.rank[x]==self.rank[y]:
#           self.rank[y] = self.rank[y]+1
   def findset(self,x):
      x=int(x)
      if x != self.p[x]:
         self.p[x]=self.findset(self.p[x])
      return self.p[x]
   def compress(self):
      nset=0
      for i in range(1,self.length+1):
         if self.p[i]==i:
            nset+=1
            self.F[i]=nset
         else:
            new=self.findset(i)
            if new<i:
               self.F[i]=self.F[new]
            else:
               print "help!",i,new,self.p[i]
@}

Now we implemement that same code in c:

Just put everything in a header file.

@O dset.h
@{
#include <stdlib.h>
#ifndef _dset_h
#define _dset_h
void dset_initialise(int size); /* array to hold real values of each */
int dset_new(void);
void dset_makeunion(int r1, int r2);
void dset_link( int r1, int r2);
int dset_find(int x);
void dset_compress(void);
int * S=NULL;
int * F=NULL;

void dset_initialise(int size){
   int i ;
   if(S==NULL){
    S=(int *) ( malloc(size*sizeof(int) ) );
    }
   else{ 
    if(S[0]!=size){S=(int *)(realloc(S,size*sizeof(int)));}
    }
   S[0]=size;
   for(i=1;i<size;i++)S[i]=0;
}

int dset_new(void){
   /* S[0] will always hold the length of the array */
   /* S[:-1] holds the current element */
   int length, current, i;
   length=S[0];
   current=(++S[S[0]-1]);
   if(current+3>length){
      S=(int *)(realloc(S,length*2*sizeof(int)));  
      S[0]=length*2;
      S[length-1]=0;
      for(i=length-1;i<length*2;i++)S[i]=0;
      S[length*2-1]=current;
   }
   S[current]=current;
   return current;
}

void dset_makeunion( int r1, int r2){
   int a,b;
   a=dset_find(r1);
   b=dset_find(r2);
   dset_link(a, b);
}

void dset_link(int r2, int r1){
   if(r1>r2){S[r1]=r2;}
   if(r1<r2){S[r2]=r1;}
   /* if r1==r2 then they are already a union */
}

int dset_find(int x){
   if (x==0){
    return 0;
   }
   if (S[x] != x){
      S[x]=dset_find(S[x]);
      }
   return S[x];
}

void dset_compress(void){
   int i,nset,new;
   if(F==NULL){
    F=(int *) ( malloc(S[0]*sizeof(int) ) );
    }
   else{
    F=(int *) ( realloc(F,S[0]*sizeof(int) ) );
    }
   for(i=0;i<S[0];i++)F[i]=0;
   nset=0;
   for(i=1;i<S[0];i++){ 
      if(S[i]==i){ 
         nset++;
         F[i]=nset;
	 }
      else{
         new=dset_find(i);
	 if (new<i){F[i]=F[new];}
	 else  {printf("i=%d new=%d S[i]=%d\n",i,new,S[i]);}
	 }
   }
}

#endif /* _dset_h */

@}


Now make a little thing to test it:

@O dset.c
@{

#include <stdio.h>
#include <stdlib.h>
#include "dset.h"
int main(int argc , char* argv[]){
    int i, j, k, Set1;
    printf("Now I am really confused");
    fflush(stdin);
    dset_initialise(100);
    printf("Got through the initialise, S=%d",S);
    fflush(stdin);
    j = k = 1;
    fflush(stdin);
    for(i=1;i<=101;i++)dset_new();
    printf("after making new, S=%d",S);    
    printf("We are set up for 200 elements");
    fflush(stdin);
    while( k <= 8 )
    {
        j = 1;
        while( j < 8 )
        {
            dset_makeunion( j, j+k );
            j += 2 * k;
        }
        k *= 2;
    }
    i = 1;
    printf("About to examine the sets");
    for( i = 1; i <= 30; i++ )
    {
        Set1 = dset_find( i );
        printf( "%d**", Set1 );
    }
    printf( "\n" );
    return 0;
}

@}

Now we have some useful little data structures for creating disjoint
sets. That is to say that if you have a pixel and want to put it in
the same set as the neighbouring pixel, this data structure will
do that.


\subsection{Collecting thresholded pixels into peaks (or objects)}

Assume we are given an image and a threshold and are asked to return
a list of peaks. The quick way to do this (apparently) is to scan 
through the image line by line. Whenever we find a pixel which is 
above the threshold we need to assign it to a peak. Checking only
the pixels which have come before (previous line and current) we
can will have some choices. If the pixel is not joined, we make a new 
peak with this pixel in it. If the pixel is joined to one other thresholded
pixel we give it the same peak label. If the pixel joins two groups
which were previously different peaks we assign the whichever we like
and make a note that those peaks are equivalent.

In python: make a new image to hold the peak assignment. Make a list of
some sort to identify the equivalences.

@D findpeaks
@{def findpeaks(image,threshold):
   """ Find pixels above threshold in image and assign them to
   peaks
   """
   peakim=zeros(image.shape,Int16)
   ds=disjointset()
   pp=0
   ppa=0
   for i in range(1,image.shape[0]-1):
      for j in range(1,image.shape[1]-1):
         if image[i,j] > threshold:
            ppa=ppa+1
            # cases are i-1,j-1/j/j+1; i,j-1
            if peakim[i-1,j-1]>0:
               peakim[i,j]=peakim[i-1,j-1]
            if peakim[i-1,j]>0:
               if peakim[i,j]==0:
                  peakim[i,j]=peakim[i-1,j]
               else:
                  ds.makeunion(peakim[i,j],peakim[i-1,j])
            if peakim[i-1,j+1]>0:
               if peakim[i,j]==0:
                  peakim[i,j]=peakim[i-1,j+1]
               else:
                  ds.makeunion(peakim[i,j],peakim[i-1,j+1])
            if peakim[i,j-1]>0:
               if peakim[i,j]==0:
                  peakim[i,j]=peakim[i,j-1]
               else:
                  ds.makeunion(peakim[i,j],peakim[i,j-1])
            if peakim[i,j]==0:
               peakim[i,j]=ds.makeset()
         if pp%1000 == 0:
            s="%d pixels processed, %d above threshold\r" % (pp,ppa)
            print s,   
         pp=pp+1
   # All pixels are now assigned to peaks, we now go through 
   # the image a second time and correct the peak indexing
   # using the disjoint set
   m=0
   ds.compress()
   for i in range(1,image.shape[0]-1):
      for j in range(1,image.shape[1]-1):
         if peakim[i,j]>0:
            v=ds.F[int(peakim[i,j])]
            peakim[i,j]=v
            if v>m:m=v
   print
   for i in range(1,m+1):
      print "i=",i,"set=",ds.F[i]
   # pixels are now assigned to peaks
   return peakim
@}

Now to test the code:

@o testpeaksearch.py
@{
import imagereaders 
@< disjointset @>
@< findpeaks @>
import blobproperties
@< testimage @>
if __name__=="__main__":
   import sys
   head,image=readimage(sys.argv[1])
#   image=testimage
   t=5000 # threshold
   from time import time
   start=time()
#   pim=findpeaks(image,t)
   end=time()
   print "\nIn pure python that took", end-start   
   cim=blobproperties.connectedpixels(image,t)
   cim2=blobproperties.connectedpixels(image[512:1024,0:512],t)
   cim3=blobproperties.connectedpixels(image[20:1004,20:1004],t)
   print "\nIn C that took", time()-end   
   print maximum.reduce(ravel(cim)),maximum.reduce(ravel(cim))
   a,x,y=blobproperties.blobproperties(image,cim)
   print "Blob area, x, y"
   for i in range(a.shape[0]):
      print i,a[i],x[i],y[i]
#   junk=raw_input("press a key to quit")
@}

Wahey. It appears the python versio worked the very first time I tried it. Staggering.
Timed at about 40 seconds for a 1K image with 4000 pixels above the
threshold, which could be significantly faster in Numeric or plain C.
Need to get the Disjoint sets thingummy implemented in C first.


\subsection{Make a test image}

For testing the blob search thing we need a test image which will give a 
unique answer for getting the algorithms right. A circle and a few 
funny shaped blobs should do it.

@d testimage
@{
testimage=zeros((512,512),Float32)
for i in range(testimage.shape[0]):
   for j in range(testimage.shape[1]):
      a=i-100
      b=j-100
      r=sqrt(a*a+b*b)
      if r > 50 and r < 75: testimage[i,j]=100
testimage[20,150:210]=100
testimage[10:200,160]=100
testimage[400:410,400:410]=100
testimage[409,350:401]=100
testimage[350:401,409]=100

@}

\subsubsection{Now to do the peak search in c}

@D jonconnect
@{
/* nx, ny are for indexing the data array only!!! */
/* if someone sends a subarray we index it like that and don't need it */
/* to be contiguous in memory */
/* dx and dy will index results array */
int i,j,k,l,m,p;
dset_initialise(16384); /* A suitably large number */
/* first row */
i=0;
/* first pixel */
if (d[0] > t) r[0]=dset_new();
else r[0]=0;
/* up to end of first row */
for(j=1 ; j<dy; j++){ 
  if (d[j*ny] > t){
     /* Assign pixel */
     if( r[j-1] > 0 ) r[j]=r[j-1];
     else r[j]=dset_new();
     }
  else r[j]=0;
  }
/* all the rest of the rows (to end) */
for(i=1 ; i<dx ; i++){
   for(j=1 ; j<dy-1; j++){    
      p=i*dx+j; /* strides for r are dx and 1, it is contiguous via malloc */
      if (d[i*nx+j*ny] > t){
         /* Assign pixel */
         k=0;l=0;
         /* i-1, j-1, assign same index */
         if( (k=r[p-dx-1]) > 0 ) {
             r[p]=k; l=k;
	     }	 
      	/* i-1,j */
         if( (k=r[p-dx]) > 0 ) {
	          if (l==0){
	             r[p]=k; l=k;
		          }
             else { dset_makeunion(k,l); }
	      }
         /* i-1, j+1 */
         if( ( k=r[p-dx+1] ) > 0 ) {
             if(l==0){
	             r[p]=k; l=k;
		          }
	          else { dset_makeunion(k,l); }
	      }
         /* i, j-1 */
         if( ( k=r[p-1] ) > 0 ) {
             if(l==0){
	             r[p]=k; l=k;
		          }
	          else { dset_makeunion(k,l); }
	          }
         if(l==0){ /* pixel has no neighbours thus far */
             r[p]=dset_new();
	     } 
	 } /* active pixel */
    else {
         r[p]=0;
	 }
  } /* j */
  j=dy-1;p=i*dx+j;
/* now do the last column with j==dy-1 */
      if (d[i*nx+j*ny] > t){
         /* Assign pixel */
         k=0;l=0;
         /* i-1, j-1, assign same index */
         if( (k=r[p-dx-1]) > 0 ) {
             r[p]=k; l=k;
	     }	 
      	/* i-1,j */
         if( (k=r[p-dx]) > 0 ) {
	          if (l==0){
	             r[p]=k; l=k;
		          }
             else { dset_makeunion(k,l); }
	      }
         /* i-1, j+1     not this one !
         if( ( k=r[p-dx+1] ) > 0 ) {
             if(l==0){
	             r[p]=k; l=k;
		          }
	          else { dset_makeunion(k,l); }
	      } */
         /* i, j-1 */
         if( ( k=r[p-1] ) > 0 ) {
             if(l==0){
	             r[p]=k; l=k;
		          }
	          else { dset_makeunion(k,l); }
	          }
         if(l==0){ /* pixel has no neighbours thus far */
             r[p]=dset_new();
	     } 
	 } /* active pixel */
    else {
         r[p]=0;
	 }

} /* i */


dset_compress();
/* Now loop through the results image again, assigning the 
                            corrected numbers to each pixel */
l=0;m=0;
for( i = 0 ; i < dx ; i++ ){    /* i,j is looping along the indices data array */
   for( j = 0 ; j < dy ; j++ ){
      p=i*dx+j;
      k=r[p];
      r[p]=F[k];
      if (r[p]>0)m++;
      if(r[p]>l)l=r[p];
      }
   }

return l ; /* number of blobs found */
/* for (i=0; i<=l; i++)printf("i=%d F[i]=%d\t",i,F[i]); */
/* printf("\nm=%d\n",m); */
}
@}

@d blobmeasure
@{
{int i, j, k, l, p;
float datum;
for(i=0 ; i<dx ; i++){
   for(j=0 ; j<dy; j++){
      p=i*dx+j;
      k=res[p]; /* strides for r are dx and 1, it is contiguous via malloc */
      /* we allow k==0 for the "not in a blob blob" */ 
      datum=dataimage[i*nx+j*ny];
      ncycles+=t2-t1;
      for(l=0;l<nprop;l++){
          blobarray[k*nprop+l]=(propertylist[l])(datum,blobarray[k*nprop+l],i,j);
         }
      }
   }
}
@}


\subsection{Accumulating mean and variance}

Thanks to the pydev mailing list summary which pointed me towards a recurrence relation
in Knuth's~\cite{Knuth::stddev} book for computing standard deviations without having 
such a bad loss
of precision. These computer scientists are frighteningly clever sometimes.

Textbooks give:

\[ \sigma = \sqrt{ \frac { 
  \left( n \sum_{1<k<n}x_{k}^{2} - \left( \sum_{1<k<n} x_{k} \right)^{2} \right) }{ n(n-1)} } \]

Apparently you can end up trying to make sqrt(negative) with this formula, essentially as
the thing is the difference of two very large numbers. The recurrence formula's are:

\[ M_{1} = x_{1} ; S_{1}=0 \]
\[ M_{k} = M_{k-1}+(x_{k}-M_{k-1})/k \]
\[ S_{k} = S_{k-1}+(x_{k}-M_{k-1})*(x_{k}-M_{k}) \]

for $ 2 \leq k \leq n$ and where $\sigma=\sqrt{S_{n}/(n-1)}$.

For the blobproperties algorithms this has some interesting niceties. Most importantly, we 
can accumulate the mean of a value directly, without doing any postprocessing after 
accumulating sums. Well, the $S_{k}$ is going to have to be divided by the number of
points.

Tim Peter's added a comment on PyDev that you do better to compute a fresh mean each
time by keeping a running sum and dividing the running sum by k. 

Sadly the variance formula needs to know the current value for the mean as well, which
contradicts the ugly way of doing this in use so far. Also I don't know how to
extend it for merging two partially completed lists. Eventually we would like to 
have this working at the disjoint set level, so as not to need to pass the data values
through the processor more than once.

I'm not actually convinced that anyone will make that much use of the variance numbers
anyway.




@O jonproperty.h
@{
#ifndef JONPROPERTY_H
#define JONPROPERTY_H
float accumulate_npt   ( float, float, int, int); /* data, current, i, j ; 1 */
float accumulate_sumI  ( float, float, int, int); /* 2 */ 
float accumulate_sumX  ( float, float, int, int); /* 3 */ 
float accumulate_sumY  ( float, float, int, int); /* 4 */ 
float accumulate_sumXI ( float, float, int, int); /* 5 */ 
float accumulate_sumYI ( float, float, int, int); /* 6 */  
float accumulate_maxI  ( float, float, int, int); /* 7 */ 
float accumulate_minI  ( float, float, int, int); /* 8 */ 
float accumulate_maxX  ( float, float, int, int); /* 9 */ 
float accumulate_minX  ( float, float, int, int); /* 10 */ 
float accumulate_maxY  ( float, float, int, int); /* 11 */ 
float accumulate_minY  ( float, float, int, int); /* 12 */ 
float accumulate_varI  ( float, float, int, int); /* 13 */ 
float accumulate_varX  ( float, float, int, int); /* 14 */ 
float accumulate_varY  ( float, float, int, int); /* 15 */ 
float accumulate_varXY ( float, float, int, int); /* 16 */ 
float accumulate_varXI ( float, float, int, int); /* 17 */ 
float accumulate_varYI ( float, float, int, int); /* 18 */  
float accumulate_varXYI( float, float, int, int); /* 19 */ 

#define NPROPERTY 19
static float (*propertylist[NPROPERTY])(float, float, int, int)={ 
accumulate_npt,
accumulate_sumI,
accumulate_sumX,
accumulate_sumY,
accumulate_sumXI,
accumulate_sumYI,
accumulate_maxI,
accumulate_minI,
accumulate_maxX,
accumulate_minX,
accumulate_maxY,
accumulate_minY,
accumulate_varI,
accumulate_varX,
accumulate_varY,
accumulate_varXY,
accumulate_varXI,
accumulate_varYI,
accumulate_varXYI } ;

static char* propertynames[NPROPERTY] ={
"accumulate_npt",
"accumulate_sumI",
"accumulate_sumX",
"accumulate_sumY",
"accumulate_sumXI",
"accumulate_sumYI",
"accumulate_maxI",
"accumulate_minI",
"accumulate_maxX",
"accumulate_minX",
"accumulate_maxY",
"accumulate_minY",
"accumulate_varI",
"accumulate_varX",
"accumulate_varY",
"accumulate_varXY",
"accumulate_varXI",
"accumulate_varYI",
"accumulate_varXYI" } ;

enum jonsproperties {npt,sumI,sumX,sumY,sumXI,sumYI,
maxI,minI,maxX,minX,maxY,minY,
varI,varX,varY,varXY,varXI,varYI,varXYI};


#endif
@}

@d jonproperty
@{
float accumulate_npt   ( float d, float c, int i, int j){ return c+1   ;}
float accumulate_sumI  ( float d, float c, int i, int j){ return c+d   ;}
float accumulate_sumX  ( float d, float c, int i, int j){ return c+i   ;}
float accumulate_sumY  ( float d, float c, int i, int j){ return c+j   ;}
float accumulate_sumXI ( float d, float c, int i, int j){ return c+i*d ;}
float accumulate_sumYI ( float d, float c, int i, int j){ return c+j*d ;} 
float accumulate_maxI  ( float d, float c, int i, int j){ if(d>c)return d;return c;}
float accumulate_minI  ( float d, float c, int i, int j){ if(d<c)return d;return c;}
float accumulate_maxX  ( float d, float c, int i, int j){ if(i>c)return i;return c;}
float accumulate_minX  ( float d, float c, int i, int j){ if(i<c)return i;return c;}
float accumulate_maxY  ( float d, float c, int i, int j){ if(j>c)return j;return c;}
float accumulate_minY  ( float d, float c, int i, int j){ if(j<c)return j;return c;}
float accumulate_varI  ( float d, float c, int i, int j){ return c+d*d; }
float accumulate_varX  ( float d, float c, int i, int j){ return c+i*i; }
float accumulate_varY  ( float d, float c, int i, int j){ return c+j*j; }
float accumulate_varXY ( float d, float c, int i, int j){ return c+i*j; }
float accumulate_varXI ( float d, float c, int i, int j){ return c+i*i*d; }
float accumulate_varYI ( float d, float c, int i, int j){ return c+j*j*d; }
float accumulate_varXYI( float d, float c, int i, int j){ return c+i*j*d; }
@}

@o jonproperty.c
@{
#include "jonproperty.h"
@<jonproperty@>
@}
@o testjonproperty.c
@{
#include <stdio.h>
#include <float.h>
#include "jonproperty.h"
#define NPTS 100

main(){

int i,j,k;
float val[NPROPERTY],d[NPTS][NPTS],r;
for(i=0;i<NPTS;i++)
   for(j=0;j<NPTS;j++){
      d[i][j]=0.0;
      r=(i-50)*(i-50)+(j-50)*(j-50)-(i-50)*(j-50);
      if(r<30)d[i][j]=500-r;
      }
for(k=0;k<NPROPERTY;k++)
   val[k]=0.;

val[minI]=FLT_MAX;
val[minX]=FLT_MAX;
val[minY]=FLT_MAX;
for(i=0;i<NPTS;i++)
   for(j=0;j<NPTS;j++)
      for(k=0;k<NPROPERTY;k++)
         if(d[i][j]>0.) val[k]=(propertylist[k])(d[i][j],val[k],i,j);
      
for(k=0;k<NPROPERTY;k++)printf("property %d named %s is %f\n",k,propertynames[k],val[k]);

printf("Mean I  %f\n",val[sumI]/val[npt]);
printf("Mean X  %f\n",val[sumX]/val[npt]);
printf("Mean Y  %f\n",val[sumY]/val[npt]);
printf("Mean XI %f\n",val[sumXI]/val[sumI]);
printf("Mean YI %f\n",val[sumYI]/val[sumI]);
printf("Var  I  %f\n",val[varI]/val[npt] -
                      (val[sumI]/val[npt])*(val[sumI]/val[npt]));
printf("Var  X  %f\n",val[varX]/val[npt] - 
                      (val[sumX]/val[npt])*(val[sumX]/val[npt]));
printf("Var  Y  %f\n",val[varY]/val[npt] - 
                      (val[sumY]/val[npt])*(val[sumY]/val[npt]));
printf("Var  XY %f\n",val[varXY]/val[npt] -
                      (val[sumX]/val[npt])*(val[sumY]/val[npt]));
printf("Var  XI %f\n",val[varXI]/val[sumI] -
                      (val[sumXI]/val[sumI])*(val[sumXI]/val[sumI]));
printf("Var  YI %f\n",val[varYI]/val[sumI] - 
                      (val[sumYI]/val[sumI])*(val[sumYI]/val[sumI]));
printf("Var XYI %f\n",val[varXYI]/val[sumI] -
                      (val[sumXI]/val[sumI])*(val[sumYI]/val[sumI]));

}
@}


This is a utility routine for getting the type of a Numeric
array for python. What I really want is a code generator
which maps numeric arrays onto c code which has the 
appropriate type definition. This has been done manually
for the peak searcher but ended up making no real difference
to processing time.

@d getval
@{
int getval(char *p, int type);

int getval(char *p, int type){
  switch (type){
     case    PyArray_CHAR   : return *(char          *)p*1;
     case    PyArray_SBYTE  : return *(signed char   *)p*1;
     case    PyArray_SHORT  : return *(short         *)p*1;
     case    PyArray_INT    : return *(int           *)p*1;
     case    PyArray_LONG   : return *(long          *)p*1;
     case    PyArray_FLOAT  : return *(float         *)p*1;
     case    PyArray_DOUBLE : return *(double        *)p*1;
#ifdef PyArray_UNSIGNED_TYPES
     case    PyArray_UBYTE  : return *(unsigned char *)p*1;
     case    PyArray_USHORT : return *(unsigned short*)p*1;
     case    PyArray_UINT   : return *(unsigned int  *)p*1;
#endif
     }
   printf("Oh bugger in getval - unrecognised numeric type\n");
   exit(1);
   return 0;
}
void ptype(int type){
  printf("Your input type was ");
  switch (type){
     case    PyArray_CHAR   : printf("PyArray_CHAR *(char *)\n");break;
     case    PyArray_SBYTE  : printf("PyArray_SBYTE *(signed char *)\n");break;
     case    PyArray_SHORT  : printf("PyArray_SHORT *(short *)\n");break;
     case    PyArray_INT    : printf("PyArray_INT *(int  *)\n");break;
     case    PyArray_LONG   : printf("PyArray_LONG *(long *)\n");break;
     case    PyArray_FLOAT  : printf("PyArray_FLOAT *(float *)\n");break;
     case    PyArray_DOUBLE : printf("PyArray_DOUBLE *(double *)\n");break;
#ifdef PyArray_UNSIGNED_TYPES
     case    PyArray_UBYTE  : printf("PyArray_UBYTE *(unsigned char *)\n");break;
     case    PyArray_USHORT : printf("PyArray_USHORT *(unsigned short*)\n");break;
     case    PyArray_UINT   : printf("PyArray_UINT *(unsigned int *)\n");break;
#endif
     }
}

@}


@O connectedpixels.c
@{
#include <Python.h>                  /* To talk to python */
#include "Numeric/arrayobject.h"     /* Access to Numeric */
#include "rdtsc.h"
/* make an image of peak assignement for pixels */

int jonconnect_c (float t, char           *d, int *r, int nx, int ny, int dx, int dy);
int jonconnect_sc(float t, signed char    *d, int *r, int nx, int ny, int dx, int dy);
int jonconnect_s (float t, short          *d, int *r, int nx, int ny, int dx, int dy);
int jonconnect_i (float t, int            *d, int *r, int nx, int ny, int dx, int dy);
int jonconnect_l (float t, long           *d, int *r, int nx, int ny, int dx, int dy);
int jonconnect_f (float t, float          *d, int *r, int nx, int ny, int dx, int dy);
int jonconnect_d (float t, double         *d, int *r, int nx, int ny, int dx, int dy);
#ifdef PyArray_UNSIGNED_TYPES
int jonconnect_uc(float t, unsigned char  *d, int *r, int nx, int ny, int dx, int dy);
int jonconnect_us(float t, unsigned short *d, int *r, int nx, int ny, int dx, int dy);
int jonconnect_ui(float t, unsigned int   *d, int *r, int nx, int ny, int dx, int dy);
#endif

void blobmeasure_c (int * res, char *             dataimage, int type, int nx, int ny, 
                                        int dx, int dy, float * blobarray, int nprop);
void blobmeasure_sc(int * res, signed char *      dataimage, int type, int nx, int ny, 
                                        int dx, int dy, float * blobarray, int nprop);
void blobmeasure_s (int * res, short *            dataimage, int type, int nx, int ny, 
                                        int dx, int dy, float * blobarray, int nprop);
void blobmeasure_i (int * res, int *              dataimage, int type, int nx, int ny, 
                                        int dx, int dy, float * blobarray, int nprop);
void blobmeasure_l (int * res, long *             dataimage, int type, int nx, int ny, 
                                        int dx, int dy, float * blobarray, int nprop);
void blobmeasure_f (int * res, float *            dataimage, int type, int nx, int ny, 
                                        int dx, int dy, float * blobarray, int nprop);
void blobmeasure_d (int * res, double *           dataimage, int type, int nx, int ny, 
                                        int dx, int dy, float * blobarray, int nprop);
#ifdef PyArray_UNSIGNED_TYPES
void blobmeasure_uc(int * res, unsigned char *    dataimage, int type, int nx, int ny, 
                                        int dx, int dy, float * blobarray, int nprop);
void blobmeasure_us(int * res, unsigned short *   dataimage, int type, int nx, int ny, 
                                        int dx, int dy, float * blobarray, int nprop);
void blobmeasure_ui(int * res, unsigned int *     dataimage, int type, int nx, int ny, 
                                        int dx, int dy, float * blobarray, int nprop);
#endif



#include "jonproperty.h"

static PyObject * getrefcount (PyObject *self, PyObject *args);

static PyObject * _connectedpixels (PyObject *self, PyObject *args, PyObject *keywds);
void initconnectedpixels(void);
#include "dset.h" 

static PyObject * getrefcount (PyObject *self, PyObject *args){
PyObject *in=NULL;
if(!PyArg_ParseTuple(args,"O",&in)) return NULL;
return Py_BuildValue("i",in->ob_refcnt-1); /* the minus 1 is for this reference here */
}


@<jonproperty@>

static PyObject * _connectedpixels (PyObject *self, PyObject *args, PyObject *keywds){
  PyArrayObject *in=NULL, *r=NULL;
  float threshold;
  int *res=NULL, size, type , nx, ny, dx, dy, i;
  int com=0, nblob=0, newdims[2];
  static char *kwlist[]={"data","threshold","Centreofmass", NULL};
  if(!PyArg_ParseTupleAndKeywords(args, keywds, "O!f|i",kwlist,   /* Pyobj, float */
                        &PyArray_Type, &in,                          /* image arg */
                        &threshold, &com))                           /* threshold */
      return NULL;
  type=in->descr->type_num;
  size=in->dimensions[0]*in->dimensions[1];
 
/* printf(
  "in->dimensions[0]=%d, in->dimensions[1]=%d, in->strides[0]=%d, in->strides[1]=%d, ",
   in->dimensions[0],in->dimensions[0],in->strides[0],in->strides[1]); */
  nx=in->strides[0];
  ny=in->strides[1];
  dx=in->dimensions[0];
  dy=in->dimensions[1];
  res =(int *) ( malloc(dx*dy*sizeof(int) ) );

  switch (type){
     case    PyArray_CHAR   : 
         /* printf("PyArray_CHAR *(char *) %d\n",sizeof(char)); */
         nx/=sizeof(char); ny/=sizeof(char);
         nblob=jonconnect_c (threshold,(char *)          in->data,res,nx,ny, dx, dy);
         break;
     case    PyArray_SBYTE  : 
         /* printf("PyArray_SBYTE *(signed char *) %d\n",sizeof(signed char)); */
         nx/=sizeof(signed char); ny/=sizeof(signed char);
         nblob=jonconnect_sc(threshold,(signed char *)   in->data,res,nx,ny, dx, dy);
         break;
     case    PyArray_SHORT  : 
         /* printf("PyArray_SHORT *(short *) %d\n",sizeof(short)); */
         nx/=sizeof(short); ny/=sizeof(short);
         nblob=jonconnect_s (threshold, (short *)        in->data,res,nx,ny, dx, dy);
         break;
     case    PyArray_INT    : 
         /* printf("PyArray_INT *(int  *) %d\n",sizeof(int)); */
         nx/=sizeof(int); ny/=sizeof(int);
         nblob=jonconnect_i (threshold, (int *)          in->data,res,nx,ny, dx, dy);
         break;
     case    PyArray_LONG   : 
         /* printf("PyArray_LONG *(long *) %d\n",sizeof(long)); */
         nx/=sizeof(long); ny/=sizeof(long);
         nblob=jonconnect_l (threshold,(long *)          in->data,res,nx,ny, dx, dy);
         break;
     case    PyArray_FLOAT  : 
         /* printf("PyArray_FLOAT *(float *) %d\n",sizeof(float)); */
         nx/=sizeof(float);ny/=sizeof(float);
         nblob=jonconnect_f (threshold, (float *)        in->data,res,nx,ny, dx, dy);
         break;
     case    PyArray_DOUBLE : 
         /* printf("PyArray_DOUBLE *(double *) %d\n",sizeof(double)); */
         nx/=sizeof(float); ny/=sizeof(float);
         nblob=jonconnect_d (threshold, (double *)       in->data,res,nx,ny, dx, dy);
         break;
#ifdef PyArray_UNSIGNED_TYPES
     case    PyArray_UBYTE  : 
         /* printf("PyArray_UBYTE *(unsigned char *) %d\n",sizeof(unsigned char)); */
         nx/=sizeof(unsigned char); ny/=sizeof(unsigned char);
         nblob=jonconnect_uc(threshold, (unsigned char *)in->data,res,nx,ny, dx, dy);
         break;
     case    PyArray_USHORT : 
         /* printf("PyArray_USHORT *(unsigned short*) %d\n",sizeof(unsigned short)); */
         nx/=sizeof(unsigned short); ny/=sizeof(unsigned short);
         nblob=jonconnect_us(threshold,(unsigned short*) in->data,res,nx,ny, dx, dy);
         break;
     case    PyArray_UINT   : 
         /* printf("PyArray_UINT *(unsigned int *) %d\n",sizeof(unsigned int)); */
         nx/=sizeof(unsigned); ny/=sizeof(unsigned);
         nblob=jonconnect_ui(threshold,(unsigned int *)  in->data,res,nx,ny, dx, dy);
         break;
#endif
     }
  if (com==0) {
     r = (PyArrayObject *) PyArray_FromDims(in->nd,in->dimensions,PyArray_INT);
     for(i=0;i<in->dimensions[0]*in->dimensions[1];i++) /* copy results array out, why? */  
        (*(int*)(r->data + i*sizeof(int) ))=res[i];
     free(res);
     return PyArray_Return(r);
  } else { /* return array of spot positions in x,y */
          /* loop over data */
     if (nblob>0){
        newdims[0]=nblob+1;
        newdims[1]=com;
        r = (PyArrayObject *) PyArray_FromDims(2,newdims,PyArray_FLOAT);
        for(i=0;i<nblob*com;i++)(*(float *)(r->data + i*sizeof(float)))=0.;
        switch (type){
        case    PyArray_CHAR   : 
         blobmeasure_c (res, (char *) in->data,type,nx,ny,dx,dy,(float *) r->data,com);
         break;
        case    PyArray_SBYTE  : 
         blobmeasure_uc(res,(signed char*)in->data,type,nx,ny,dx,dy,(float*)r->data,com);
         break;
        case    PyArray_SHORT  : 
         blobmeasure_s(res, (short *) in->data,type,nx,ny,dx,dy,(float *) r->data,com);
         break;
        case    PyArray_INT    : 
         blobmeasure_i(res, (int *) in->data,type,nx,ny,dx,dy,(float *) r->data,com);
         break;
        case    PyArray_LONG   : 
         blobmeasure_l(res, (long *) in->data,type,nx,ny,dx,dy,(float *) r->data,com);
         break;
        case    PyArray_FLOAT  : 
         blobmeasure_f(res, (float *)in->data,type,nx,ny,dx,dy,(float *) r->data,com);
         break;
        case    PyArray_DOUBLE : 
         blobmeasure_d(res, (double *) in->data,type,nx,ny,dx,dy,(float *) r->data,com);
         break;
#ifdef PyArray_UNSIGNED_TYPES
       case    PyArray_UBYTE  : 
        blobmeasure_uc(res,(unsigned char*)in->data,type,nx,ny,dx,dy,(float*)r->data,com);
        break;
       case    PyArray_USHORT : 
        blobmeasure_us(res,(unsigned short*)in->data,type,nx,ny,dx,dy,(float*)r->data,com);
        break;
       case    PyArray_UINT   : 
        blobmeasure_ui(res,(unsigned int*)in->data,type,nx,ny,dx,dy,(float*)r->data,com);
        break;
#endif
     }
        free(res);
        return PyArray_Return(r);
    }
    else { free(res);return NULL; }
  }


}



void blobmeasure_c (int * res, char *             dataimage, int type, int nx, int ny, 
                                        int dx, int dy, float * blobarray, int nprop)
@< blobmeasure @>
void blobmeasure_sc(int * res, signed char *      dataimage, int type, int nx, int ny, 
                                        int dx, int dy, float * blobarray, int nprop)
@< blobmeasure @>
void blobmeasure_s (int * res, short *            dataimage, int type, int nx, int ny, 
                                        int dx, int dy, float * blobarray, int nprop)
@< blobmeasure @>
void blobmeasure_i (int * res, int *              dataimage, int type, int nx, int ny, 
                                        int dx, int dy, float * blobarray, int nprop)
@< blobmeasure @>
void blobmeasure_l (int * res, long *             dataimage, int type, int nx, int ny, 
                                        int dx, int dy, float * blobarray, int nprop)
@< blobmeasure @>
void blobmeasure_f (int * res, float *            dataimage, int type, int nx, int ny, 
                                        int dx, int dy, float * blobarray, int nprop)
@< blobmeasure @>
void blobmeasure_d (int * res, double *           dataimage, int type, int nx, int ny, 
                                        int dx, int dy, float * blobarray, int nprop)
@< blobmeasure @>
#ifdef PyArray_UNSIGNED_TYPES
void blobmeasure_uc(int * res, unsigned char *    dataimage, int type, int nx, int ny, 
                                        int dx, int dy, float * blobarray, int nprop)
@< blobmeasure @>
void blobmeasure_us(int * res, unsigned short *   dataimage, int type, int nx, int ny, 
                                        int dx, int dy, float * blobarray, int nprop)
@< blobmeasure @>
void blobmeasure_ui(int * res, unsigned int *     dataimage, int type, int nx, int ny, 
                                        int dx, int dy, float * blobarray, int nprop)
@< blobmeasure @>
#endif

int jonconnect_c (float t, char           *d, int *r, int nx, int ny, int dx, int dy){
@< jonconnect @>
int jonconnect_sc(float t, signed char    *d, int *r, int nx, int ny, int dx, int dy){
@< jonconnect @>
int jonconnect_s (float t, short          *d, int *r, int nx, int ny, int dx, int dy){
@< jonconnect @>
int jonconnect_i (float t, int            *d, int *r, int nx, int ny, int dx, int dy){
@< jonconnect @>
int jonconnect_l (float t, long           *d, int *r, int nx, int ny, int dx, int dy){
@< jonconnect @>
int jonconnect_f (float t, float          *d, int *r, int nx, int ny, int dx, int dy){
@< jonconnect @>
int jonconnect_d (float t, double         *d, int *r, int nx, int ny, int dx, int dy){
@< jonconnect @>
#ifdef PyArray_UNSIGNED_TYPES
int jonconnect_uc(float t, unsigned char  *d, int *r, int nx, int ny, int dx, int dy){
@< jonconnect @>
int jonconnect_us(float t, unsigned short *d, int *r, int nx, int ny, int dx, int dy){
@< jonconnect @>
int jonconnect_ui(float t, unsigned int   *d, int *r, int nx, int ny, int dx, int dy){
@< jonconnect @>
#endif




static PyMethodDef _connectedpixelsMethods[] = {
   {"_connectedpixels", (PyCFunction) _connectedpixels, 
    METH_VARARGS | METH_KEYWORDS ,
   "Assign connected pixels in image (first arg) above threshold (second arg)"},
   {"getrefcount", getrefcount, METH_VARARGS, "Return the objects reference count"},
   {NULL, NULL, 0, NULL} /* setinel */
};

void init_connectedpixels(void)
{
   PyObject *m, *d;

   m=Py_InitModule("_connectedpixels", _connectedpixelsMethods);
   import_array();
   d=PyModule_GetDict(m);
}
@}


@o testconnectedpixels.py
@{
import imagereaders
import _connectedpixels
if __name__=="__main__":
   import sys
   head,image=readid11edf(sys.argv[1])
   from time import time
   start=time()
   peakim=_connectedpixels._connectedpixels(image.astype(Float32),3000)
   end=time()
   print "\nThat took", end-start
   #from imview import view
   #view(peakim,(512,512))
   peakim=_connectedpixels._connectedpixels(image.astype(Float32),3000,Centreofmass=6)
   print peakim.shape, "that took", time()-end
   print peakim[0:10,:]
#   junk=raw_input("press a key to quit")
@}



\chapter{Crystal structure refinement}

Initially take the cctbx routines and decide later if it 
is worth replacing them with anything else.
We will also have globular type objects, which need adding
onto structure factors as well, so we can decompose a 
structure factor into the sum of a normal atomic contribution
and a part due to globular objects.

\section{Globular objects}


\subsection{A hard sphere}
From a book on polymers and neutron scattering which I really
should have noted a reference for... 
The form factor for a sphere is given by the fourier
transform of the scattering density.
So that if we normalise to a unit volume we get:
\[ f(q) = \frac{1}{V} \int_V e^{-i\mb{qr}} dV = 
          \frac{1}{V} \int_V e^{-i\mb{qr}} r^2 \sin{\theta} d\theta d\phi dr \]
where the volume element $dV =  r^2 \sin{\theta} d\theta d\phi dr $
is given by $dr$ as a shell radius, $rd\theta$ as the width over
the $\theta$ in radians angle (range $0$ to $\pi$) and
$ r\sin{\theta} d\phi$ for the rotation in the $xy$ plane, if
$\theta=0$ is along $z$.
In the complex exponential we have the dot product $\mb{qr}$, which
ought to be given by $|r||q|\cos{\theta}$ if $\mb{q}$ is along
the $z$ axis.
For convenience a new variable $u=\cos\theta$ is used and the 
integral becomes:
\[ f(q) = \frac{1}{V} \int_{r=0}^{r=R} \int_{\phi=0}^{\phi=2\pi}
\int_{u=-1}^{u=1} e^{-iqru} r^2 dr d\phi du \]
\[ f(q) = \frac{2\pi}{V} \int_{r=0}^{r=R} 
\int_{u=-1}^{u=1} e^{-iqru} r^2 dr du \]
Now we need to integrate over u to recover a ``well known result''. 
Hilarious... 
\[ f(q) = \frac{2\pi}{V} \int_{r=0}^{r=R} \left[\frac{e^{-iqru}}{iqr} 
      \right]_{u=-1}^{u=1} r^2 dr \]
\[ f(q) = \frac{2\pi}{V} \int_{r=0}^{r=R} \left[
     \frac{e^{-iqr} - e^{iqr}}{iqr} \right] r^2 dr \]
\[ f(q) = \frac{2\pi}{V} \int_{r=0}^{r=R} \left[
     \frac{\cos(-qr)-\cos(qr)+i(\sin(-qr)-\sin(qr)}{iqr} \right] dr \]
Now $\sin(-x)=-\sin(x)$ and $\cos(-x)=cos(x)$ as they are odd or
even function respectively, so that we get
\[ f(q) = \frac{2\pi}{V} \int_{r=0}^{r=R}\frac{2\sin(qr)}{qr}  r^2 dr \]
The volume of a sphere is given by $4\pi R^3/3$ so we can get rid
of the two factors or two and replace the $V$ and $\pi$ by
$R^3$.
\[ f(q) = \frac{4\pi}{V=\frac{4\pi R^3}{3}}\int_{r=0}^{r=R}\frac{\sin(qr)}{qr}
   r^2 dr \]
\[ f(q) = \frac{3}{R^3}\int_{r=0}^{r=R}\frac{\sin(qr)}{qr} r^2 dr \]
\[ f(q) = \frac{3}{qR^3} \int_{r=0}^{r=R}\sin(qr) r dr \]
That was the ``well known result''. Love a sense of humour.
Now the integration over $r$ is done by parts, which means that:
\[ \int_a^b u dv = \left| uv \right|_a^b - \int_a^b v du \]
which is something like the integral way of using the product rule
for differentiating. 
We are hoping to split up the integral into two things, one of which
is something we can integrate and one which is something we can 
differentiate.
If we make $u=r$ and $dv=\sin(qr) dr$ then we will need
$du=dr$ and $v=\cos(qr)/q$.
\[ f(q) = \frac{3}{qR^3} \left[ uv - \int v du \right]  \]
\[ f(q) = \frac{3}{qR^3} \left[ 
    \left| \frac{r \cos(qr)}{q} \right|_{r=0}^{r=R} - % uv
     \int_{r=0}^{r=R}  \frac{\cos(qr)}{q} dr \right] \]
\[ f(q) = \frac{3}{qR^3} \left[ 
   \frac{\sin(qR)}{q^2}  - \frac{qR \cos(qR)}{q^2} \right] \]
\[ f(q) = \frac{3}{(qR)^3} \left( \sin(qR) - qR \cos(qR) \right) \]
So, finally, you have it. 
The result that most people just directly quote. School was  a
long time ago now!

\subsection{Shells}

Not seashells, but radial shells of scattering density.
If the sphere does not have a uniform density, but instead
it has some density which is a function of radius (only) then
we can try to generalise one of the equations above.
Since we still have spherical symmetry, it makes most sense to 
skip down to the part where we integrated over $r$. 
The equation in question is:
\[ f(q) = C \int \frac{ n(r) r \sin(qr)}{q} dr \]
where $C$ is a normalisation constant defined such that
\[ \frac{1}{C} = \int n(r) r^2 dr \]
making $f(0)=1$.
For a shell of internal radius $R_i$ and external radius $R_e$
we are going to have the same integral as before, but going
from $R_i$ to $R_e$ instead of $0$ to $R$.
\[ f(q) = \frac{3}{q^3} \left[
\frac{\sin(qR_i) - qR_i\cos(qR_i)}{R_i^3}
- \frac{\sin(qR_e) - qR_e\cos(qR_e)}{R_e^3}
\right] \]
If we make the thickness of the shell small then $R_i \approx R_e \approx R$
and I've no idea how to look for the limiting behaviour of that equation.
However, putting a delta function into the integral above would seem 
to say we can just plop out the integral sign giving:
\[ f(q) = C \frac{sin(qR)}{q} \]
where the factor $C$ will tidy up any details for us. 

What would be interesting here is to consider various functional 
forms for the $n(r)$ bit to see if we can look at density as
a function of distance.
Means we need to think of things which are integrable analytically
when multiplied by $r\sin(qr)$ and which look like some kind of a 
peak. 
We could have a multishell model, where the density at each radius is
given by a constant number. 
How would we go about determining the numbers for each shell??

\subsection{Other shapes}

Eventually cylinders/rods, rings and any other blobs.
Try to make a continuous description which goes from 
isolated atoms to a delocalised density. 
That means writing the localised atoms structures in 
terms of spherical harmonic thingy.




\section{The structure factor}

Outline the mathematics and derive the derivatives of $F(hkl)$
with respect to $x,y,z,f',f''$ and solvent scattering.
Indicate how a derivative of $F(hkl)$ becomes a derivative of 
$I(hkl)$.

The structure factor is given by:
\[ F(hkl) = \sum_{i} f_i \exp(2\pi i(hx_i + ky_i + lz_i))T \]
Derivatives with respect to x,y,z:
\[ \frac{\partial F(hkl)}{\partial x} = 2\pi i h f_i \exp(2\pi i(hx + ky + lz))T \]
\[ \frac{\partial F(hkl)}{\partial x} = 2\pi h f_i \left[
                  \cos(2\pi(hx + ky + lz)) +
              i   \sin(2\pi(hx + ky + lz)) \right] T \]
If $f_i=f_0 + f' + if''$ then:
\[ \frac{\partial F(hkl)}{\partial f'} =    \exp(2\pi i(hx + ky + lz))T \]
\[ \frac{\partial F(hkl)}{\partial f''} = i \exp(2\pi i(hx + ky + lz))T \]
If we define solvent scattering in terms of the same model as is present 
in GSAS then:
\[ f_i=f_0 + f' + if'' - A \exp(-B q^2)T \]
\[ \frac{\partial F(hkl)}{\partial A} = - \exp(-B q^2)  \exp(2\pi i(hx + ky + lz))T \]
\[ \frac{\partial F(hkl)}{\partial B} = - AB \exp(-B q^2)  \exp(2\pi i(hx + ky + lz))T \]
The temperature factor has just been left in as a constant $T$ above.

These are one atom term derivatives.
If atoms are related by a symmetry operator then we need to generate
the new atomic position using the symmetry operator and convert
the derivatives with respect to x,y,z back again via the
inverse of the symmetry operator (I think?).

\section{The thermal factor}

List the various unit used for anisotropic (and isotropic) thermal
factors.
Derive the derivatives of $F(hkl)$ with respect to $U_{ij}$ and determine
how they add up when the atom is transformed to symmetry equivalent
positions.
This overlaps with the symmetry constraints below a little bit?

The temperature factor can be given by:
\[ T(hkl) = \exp{ - \mb{H B H} }\]
where B is a matrix.
So multiplying out you get:
\[ \mb{HBH} = B_{11}h^2 + B_{22}k^2 + B_{33}l^2 + B_{12}hk + B_{13}hl + B_{23}kl \]
In order to have the thermal factors in real space angstroms most authors
normalise the $\mb{H}$ vector to be in terms of a reciprocal space distance
in reciprocal angstroms ($U_{ij}$).

Derivatives are therefore:
\[ \frac{\partial T(hkl)}{\partial B_{ij}} = -\mb{H}_i \mb{H}_j T(hkl) \]
with whatever multiplier you need to convert the units of $B_{ij}$ into
the ones you are using.

\section{Structure factors from cctbx}

While it might eventually be nice to be computing structure factors
and derivatives in code here, for now we can leverage off of the
work done by others and just use the routines from cctbx.

@o cctbxpdbsf.py
@{
from cctbx import xray, crystal, miller, eltbx
from iotbx.pdb.xray_structure import from_pdb

import sys

structure = from_pdb(sys.argv[1])

structure.show_summary()

f_calc = structure.structure_factors(d_min=10.,)
fc=f_calc.f_calc()
fc.show_summary().show_array()
@}


If we want to use this without the trauma of learning how to use the cctbx array
functionality a little bit of code might come in handy:


@o cctbxfcalc.py
@{
from cctbx import xray, crystal, miller, eltbx
from cctbx.array_family import flex
from cctbx.xray import ext
import cctbx
import Numeric
def getfcalc(peaks,s,cs=None):
   """
   peaks = hkl peak list
   s  = structure
   cs = crystal symmetry
   """
   if cs==None:
      sym=s.space_group().type().lookup_symbol()
      cel=s.unit_cell().parameters()
      cs=crystal.symmetry(unit_cell=cel,space_group_symbol=sym)
   mi=flex.miller_index(peaks)
   data=flex.double((0,)*mi.size())
   ms=miller.set(cs,mi)
   ma=miller.array(ms,data)
   fcm=xray.structure_factors.from_scatterers(crystal_symmetry=cs,d_min=0.1)
   fc=fcm(s,ma,"direct").f_calc()
   f=flex.to_list(fc.data())
   return Numeric.array(f)

def getfcalcfrompdb(pdbfile,peaks=None):
   from iotbx.pdb.xray_structure import from_pdb
   s = from_pdb(pdbfile)
   if peaks==None:
       # generate hkl peaks from cctbx?
       raise Exception("Need to generate default peaks")
   return getfcalc(peaks,s)

@}




\section{Rigid molecular objects}

For doing molecular replacement it is useful to be able to make a rigid
body refinement of a protein without always expressing the structure
factor in terms of the consituent atoms.
By examining some structure factor algebra you will see that moving
a molecule by some displacement adds only a phase factor to the 
F(hkl) values. 
When there are several symmetry related molecules, the phase factor
will generally be different for each molecule.
So translational refinements can be done by making a set of structure
factors for each molecule (in space group P1), multiplying them by the
appropriate phase factor and then adding them together.
This can be checked by making the equivalent translation to a single 
molecule and putting it in the correct space group.

Rotations are generally more difficult - it involves moving the lattice 
points in reciprocal space. 
If we can assume that the molecular transform is a well behaved function
we can compute the structure factors for a single molecule in a relatively
large unit cell (eg: several times the size of the real unit cell).
Then plot these hkl values in reciprocal space.
When we want to find the hkl values for a particular orientation of the 
molecule we need to interpolate the hkl values computed for the big unit
cell. We will also need to find out how to recover the rotation and
translation operations which are being carried out in reciprocal space
to return to real space.

\subsection{Translational components}

Consider a molecular entity with fixed orientation having co-ordinates
$(x,y,z)$, it has structure factor:
\[ F_0(hkl) = \sum_{atoms} f_i \exp{2\pi i(hx+ky+lz)} \]
If it is translated to $(x+dx,y+dy,z+dz)$ the structure factors become:
\[ F_x(hkl) = \sum_{atoms} f_i \exp{2\pi i[h(x+dx)+k(y+dy)+l(z+dz)]} \]
\[ F_x(hkl) = \sum_{atoms} f_i \exp{2\pi i (hx+ky+lz)} 
                             \exp{2\pi i(h.dx+k.dy+l.dz)} \]
For all atoms the term $h.dx+k.dy+l.dz$ is equal, so that it can be 
factored out of the summation over atoms and we get:
\[ F_x(hkl) = F_0(hkl) \exp{ 2 \pi i (h.dx + k.dy + l.dz) } \]
This works for one single molecule in the unit cell.
When the unit cell contains more than one molecule we need
to compute the $F_0(hkl)$ values for each of the symmetry related
molecules and apply the symmetry operation to the vector $dx,dy,dz$ (only
the multiplicative part).

So if we have a space group which is defined by a list of symmetry operators
we compute the structure factor by making a P1 structure factor for
each of the molecules generated by symmetry. 
Then add the translational component adjusted for symmetry to each molecule.
Then add together the contributions from the different molecules.

In order to optimise the position of the molecule we will need derivatives
of the structure factors with respect to the position of the molecule.
\[ \frac{dF_x(hkl)}{d (dx)} =
   F_0(hkl) \frac{d}{d (dx)}  \exp{ 2 \pi i (h.dx + k.dy + l.dz) } \]
\[ \frac{dF_x(hkl)}{d (dx)} =
  2 \pi i h F_0(hkl)   \exp{ 2 \pi i (h.dx + k.dy + l.dz) } =
  2 \pi i h F_x(hkl)  \]
\[ \frac{dF_x(hkl)}{d (dy)} =  2 \pi i k F_x(hkl)  \]
\[ \frac{dF_x(hkl)}{d (dz)} =  2 \pi i l F_x(hkl)  \]

So we should be ready to put something together for rigid body optimisation
with respect to position of the molecule now.

\subsection{Orientational components}

Assuming we have put the molecule into spacegroup P1 and computed
the structure factors on an effectively infinite grid then we
should have the molecular transform.
To find the structure factors for some specific orientation we then
need to rotate the molecular transform to our new orientation
and lookup the values.
Which co-ordinate system is used for the transformation might
eventually be important, but for now we will just note that
all rotations can be represented by a rotation matrix
and we will expect a suitable matrix to be specified.


\section{Solvent scattering}

There are a variety of more or less complex models for the 
solvent scattering which comes for the large "void" areas
in macromolecular compounds.
We begin with the least sophisticated:

\subsection{Simplest model}

We assume the solvent scattering is in opposition to the molecular
scattering (it is not) and modify the structure factors by a 
smooth curve.
This zeroth order model is expected to fail but provide a 
reasonable first order approximation....

\[ F_{scalc} = F_{calc} [ 1 - A e^{-B s^2}] \]

where $s=\sin(\theta/\lambda)$. $A$ represents the electron
density of the solvent and $B$ represents a "thermal" factor
for it - as it mainly contributes at low angles.

If we have a list of peaks we need to compute the d-spacing
or $\sin(\theta/\lambda)$ for each one. Then

\[ \frac{d F_{scalc}}{d A} = - F_{calc} e^{-B s^2} \]
\[ \frac{d F_{scalc}}{d B} =  s^2 A F_{calc} e^{-B s^2} \]

Or for intensity:

\[ I = |F|^2 = |F_{scalc}|^2 [1 - A e^{-Bs^2}]^2  = |F_{calc}|^2 [1 - 2 A e^{-Bs^2} + A^2 e^{-2Bs^2}] \]
\[ \frac{dI}{dA} =  |F_{calc}|^2  [- 2 e^{-Bs^2} + 2 A e^{-2Bs^2}] \]
\[ \frac{dI}{dB} =  |F_{calc}|^2  [2As^2e^{-Bs^2} - 2s^2A^2e^{-2Bs^2}   ] \]

...and we could check that by using the chain rule:

\[ \frac{dI}{dA} = \frac{dI}{dF}\frac{dF_{scalc}}{dA} = 2 F_{scalc}[ - F_{calc} e^{-B s^2}] \]
\[ \frac{dI}{dA} = 2 F_{calc} [ 1 - A e^{-B s^2} ][ - F_{calc} e^{-B s^2}] \]
\[ \frac{dI}{dA} = - 2 |F_{calc}|^2 [{ e^{-Bs^2} - A e^{-2 B s^2} }] \]
\[ \frac{dI}{dB} = \frac{dI}{dF}\frac{dF_{scalc}}{dB} = 2 F_{scalc}[ s^2 A F_{calc} e^{-B s^2}  ] \]
\[ \frac{dI}{dB} = 2 F_{calc} [ 1 - A e^{-B s^2} ]   [ s^2 A F_{calc} e^{-B s^2}  ] \]
\[ \frac{dI}{dB} = 2 |F_{calc}|^2 [ s^2 A e^{-B s^2} - s^2 A^2 e^{-2B s^2} ] \]

So now if we have a list of structure factors for a particular computed model
we should be able to attempt to refine an overall scale factor and a solvent
scattering model.

\subsection{Test refinement of solvent scattering}

If we have a dils file from mprodd and a crystal structure we can attempt 
to optimise the solvent scattering parameters. This is the first useful
thing to be done with linarp during the historical course of writing
it. 

We will make a little command line script which takes the name of a dils
file, the name of a pdb file, and initial values for the solvent scattering.
It will try to refine the overall scale and two solvent scattering parameters.

@o solventmodel.py
@{
@< pycopyright @>
from Numeric import *
import model
class solventmodel(model.model):
   def __init__(self,**kwds):
      self.variables=['scale','Asolv','Bsolv']
      model.model.__init__(self,**kwds)
   def compute(self,data):
      # peak here needs to specify h,k,l,fcalc
      self.ycalc=zeros(data.npoints,Float32)
      self.gradients=zeros((data.npoints,3),Float32)
      scale=self.vv[self.vd['scale']]
      A =   self.vv[self.vd['Asolv']]
      B =   self.vv[self.vd['Bsolv']]
      for i in range(data.npoints):
         s = data.sinthetaoverlambda2(i)   # peaks must provide their d-spacings and icalc
         fsq=data.fsq[i]
         self.ycalc[i]=scale * fsq * ( 1.0 - 2.* A * exp(-B * s ) + 
                                             A * A * exp(-2.* B * s ) ) 
         self.gradients[i,0]= self.ycalc[i]/scale  # Scale factor
         self.gradients[i,1]= 2 * scale * fsq * (A*exp(-2*B*s)-exp(-B*s))
         self.gradients[i,2]= 2 * scale * s * A * fsq * (exp(-B*s)-A*exp(-2*B*s))
#         print data.ciiobject.HKL[i],s,self.ycalc[i],self.gradients[i,:]
   def gradient(self, variable):
      if variable=='scale':
         return self.gradients[:,0]
      if variable=='Asolv':
         return self.gradients[:,1]
      if variable=='Bsolv':
         return self.gradients[:,2]

if __name__=="__main__":
   # test - read the command line arguments and use the functions
   import sys
   import solventrefinedata
   import densemodelfit
   try:
      ciifile = sys.argv[1]
      cclfile = sys.argv[2]
      pdbfile = sys.argv[3]
      Asolv   = float(sys.argv[4])
      Bsolv   = float(sys.argv[5])

   except:
      print "Usage: %s ciifile cclfile pdbfile A_solv B_solv "%(sys.argv[0])
      print "Attempts to fit solvent scattering parameters"
      sys.exit()
   
   ciidataobj=solventrefinedata.solventrefinedata(ciifile,cclfile,pdbfile)
   # estimate the scale factor:
   calc=sum(ciidataobj.fsq)
   obs= sum(ciidataobj.ciiobject.Ihkl)
   s = obs/calc
   model = solventmodel(Asolv=Asolv,Bsolv=Bsolv,scale=s)

   optimiser = densemodelfit.densemodelfit(model=model,data=ciidataobj,marq=2.0)
   
   optimiser.refine()
   while raw_input("More cycles?") not in ["Q","q","n","N","0",""]:
      optimiser.refine()
@}


\section{Structure refinement}

We will begin with a very simple example model, which has already been successfully
refined using correlated integrated intensities in order to check the 
implementation of the optimiser and maths.
Later we will need to make the crystal structure model object somewhat more
efficient in terms of applying shifts etc.

FIXME: Put the LaB6 model in here.

\section{Magnetic scattering}

No immediate plans to include magnetism in the program, but eventually
it has to be there.
Whatever is done has to be sufficiently generalised to make it worth
doing magnetism at all - which means using the propagation vector
to index any magnetic peaks.

Need to be able to refine the propagation vector (gives derivatives on
peak positions).
Any implementation of magnetism \emph{must} use the symmetry of the
crystal structure. 
None of this expanding to P1 crap - either you do the symmetry properly
or you don't make it part of this program please!

\chapter{Symmetry constraints}

Derive the contraints on parameters determined by symmetry. 
These include, but should not be limited to:
\begin{itemize}
\item Unit cell parameters
\item Peak width parameters
\item hkl systematic absences
\item hkl peaks equivalent to each other
\item Factor for making normalised structure factors ($E(hkl)$)
\item Any phase relationships? eg - equal intensity, different phases?
\item Atomic positions
\item Atomic thermal factors (polar vectors)
\item Magnetic moments (axial vectors)
\end{itemize}

Read up on what is available in cctbx for carrying all of this out.
I think Ralf Grosse-Kunstleve has implemented most of the things
you can dream of doing with symmetry already?

\chapter{Restraints}

Lots of things need to be restrained for protein refinements.
The MMTK package for python implements an amber force field.
It should be possible to use that package for doing a first
attempt at restraints. 
Just need to figure out either how to supply derivatives of the 
$\chi^2_I$ function to MMTK or the other way around and get 
derivatives with respect to co-ordinates.
Also interesting is to derive normal modes for the protein
and try refinement as rigid body, followed by normal modes,
followed by restrained co-ordinates.

Generally we need the following things to be restrainable:
\begin{itemize}
\item Bond distances (2 atoms)
\item Bond angles (3 atoms)
\item Torsion angles (4 atoms)
\item Planarity (at least 4 atoms, unlimited in general)
\item Ramachandran - these are torsion angle restraints?
\item Chemical composition (constrain also)
\item Anisotropic thermal factors (TLS, constraint option)
\end{itemize}

\chapter{Corrections}

Allow the various things to be chained together to
make the final computed thing by putting them in
some sort of list.
\[I_{c,\mb{H}} = |F_{H}|^{2} A_{\mb{H}} B_{\mb{H}} C_{\mb{H}} ... \]
\[\frac{\partial I_{c,\mb{H}} }{\partial A} = \frac{I_{c,\mb{H}}}{A} \]
\[\frac{\partial I_{c,\mb{H}} }{\partial x_A} =
 \frac{\partial I_{c,\mb{H}} }{\partial A} \frac{\partial A}{\partial x_A} \]
So if there can be an arbitrarily long list of things being 
multiplied together to make a product. 
The derivative of the final product with respect to something
which varies one of the items in the list is the final product divided
by the item, multiplied by the derivative in of the item in the list.
The equations are hopefully more clear.

FIXME: Check if this is how CALPRD in PRODD actually does it.


\section{Lorentz and Polarisation}

Look up the various lorentz and polarisation factors for powder 
diffraction. 
The polarisation is a refinable parameter, and can also be expressed
in terms of the monochromator angle for laboratory diffractometers.

\section{Extinction}

Look up the Sabine paper about extinction in powders.
The refineable parameter is crystallite size, but I think
there were equations for Lorentzian and Gaussian something or other.
Check if the extinction should vary across the width of a peak - 
is it coupled with peakshape?

\section{Absorbtion}

Various possibilities, some with refineable parameters, some without.

\begin{itemize}
\item Allow the user to read in a tabulated absorbtion correction.
\item Cylindrical sample, debye-scherrer, small $\mu R$
\item Cylindrical sample, debye-scherrer, large $\mu R$
\item Thin walled cylinders
\item Flat plate
\item Flat plate with surface roughness
\end{itemize}

\section{Preferred orientation}

Need to know the scattering vector associated with each hkl
reflection in order to be able to compute the preferred orientation
corrections.
For some types of data this can be computed from the d-spacing.
For proper texture experiments the data will actually have to
provide more than just the hkl of the peaks.

Various levels of approximation are possible.
March-Dollase is something like an ellipsoid in reciprocal space.
Could think about trying to use the same functional descriptions
for PO as for peak shape effects (higher and higher rank tensors).

Generation of pole figures?

Spherical harmonics?



\chapter{Derived parameters and error estimates}

One of the major goals for this project is to learn how
to get error bars correct. 
There are two aspects to doing this - the first and most difficult
is learning how to propagate a variance-covariance matrix
from one set of parameters to another, via a non-linear function.
I would like to implement this in a fairly general way, with the hope
that it can get the answers right (discovering for itself that the
esd on $\gamma=120^\circ$ is zero for a hexagonal crystal).
The second and more interesting is to work out what the variance-
covariance matrix should actually be in the cases where the 
model does not fit the data.

\section{Error propagation}

From ``mathworld'' there are some notes on Error Propagation, which 
I am copying out here in order to try to get it straight in my head.
So:
Given $y=f(x)$ with an error of $dx$ in $x$, the absolute error
in $y$ is $dy$. The relative error is $dy/y$. 
If $x=f(u,v)$ then 
\[ x_i - <x> = (u_i -<u>) \frac{\partial x}{\partial u} +
               (v_i -<v>) \frac{\partial x}{\partial v} +...,\]
where $<x>$ denotes the mean so that:
\begin{eqnarray}
 \sigma_x^2 & = &  \frac{1}{N-1}\sum_{i=1}^N (x_i - <x>)^2 \\
            & = & \frac{1}{N-1}\sum_{i=1}^N \left[
              (u_i - <u>)^2 \left( \frac{\partial x}{\partial u}\right)^2 
              +(v_i - <v>)^2 \left( \frac{\partial x}{\partial v}\right)^2 
                \right. \\
      &   & \left. 
            +2(u_i-<u>)(v_i-<v>)\left( \frac{\partial x}{\partial u}\right)
                                \left( \frac{\partial x}{\partial v}\right)
            + ... \right] \\
\end{eqnarray}
Variance and covariance are then defined as:
\[ \sigma^2_u = \frac{1}{N-1} \sum_{i=1}^{N} (u_i-<u>)^2 \]
\[ \sigma^2_v = \frac{1}{N-1} \sum_{i=1}^{N} (v_i-<v>)^2 \]
\[ \sigma^2_{uv} = \frac{1}{N-1} \sum_{i=1}^{N} (u_i-<u>)(v_i-<v>) \]
where $\sigma_{ii}=\sigma_i^2$ so that:
\[ \sigma_x^2 = \sigma_u^2\left( \frac{\partial x}{\partial u}\right)^2
 + \sigma_v^2   \left( \frac{\partial x}{\partial v}\right)^2
 + 2\sigma_{uv} \left(\frac{\partial x}{\partial u}\right)
                \left(\frac{\partial x}{\partial v}\right)
 + ...\]

The citation for that mathematics is: ``Eric W. Weisstein. Error Propagation.
From Mathworld http://mathworld.wolfram.com/ErrorPropagation.html.
They run through a few examples using logarithms, exponentials, sum 
multiply and divide operations.

Need to have a bit more of a think about this.





\section{Error estimation}

For estimating errors when the model fits the data, we just
need to transform the error bars on the data into error
bars on the model via an appropriate transformation matrix
operation.
When the model does not fit the data we have the problem that
the transformation of ones set of error bars into another
makes no sense at all - the observations don't agree with
the model.
One approach to sorting out the problem is to increase the
error bars on the data until the model does fit. 
In many fitting programs this is implemented by multiplying
the matrix by a $\chi^2$ factor such that the $\chi^2$ for
your fit would come out to be 1.

Sadly the $\chi^2$ factor makes little sense for certain
cases which can frequently arise in powder diffraction.
For example - you might be fitting a pattern which contains
many background points and a few narrow peaks.
If the crystal structure is the problem then this only affects
the points which have peaks contributing to them. 
This leads to the general idea that it is the points which 
are wrong which are the ones that need to be corrected.
To make the assumption about transforming error bars from
data to model, we should find that the distribution of weighted
differences between model and data is also Gaussian.
That is to say that we would adjust the weights on the datapoints
such that a normal distribution is followed.

\subsection{Generating a normal probability plot}

We need to predict the weighted differences for a particular number of
uncorrelated datapoints.
If we plot a histogram of the observed weighted differences this 
should be something which can be fitted to a Gaussian.
If we plot a Gaussian and divide the area under it into the number 
of points observed, then we have the centre of each bin as the expected
weighted difference for that datapoint.

A normalised Gaussian is given by:
\[ P(x,\mu,\sigma) = \frac{1}{\sqrt{2\pi}\sigma} 
                    \exp{ \left(\frac{x-\mu}{4\sigma} \right)^2}\]
We need the integral underneath a Gaussian to get the probabilities,
and unfortunately there is no closed form analytical solution
to this problem.
Happily a fortran routine for calculating the error function
exists, so we will just need to use it to figure out where
the boundaries should be placed for the expected datapoints.
When there are a large number of datapoints it will make more
sense to do this in fortran as well.

@O normp.f
@{
c Found on the web at:
c http://lib.stat.cmu.edu/apstat/66
c
        SUBROUTINE NORMP(Z, P, Q, PDF)
C
C       Normal distribution probabilities accurate to 1.e-15.
C       Z = no. of standard deviations from the mean.
C       P, Q = probabilities to the left & right of Z.   P + Q = 1.
C       PDF = the probability density.
C
C       Based upon algorithm 5666 for the error function, from:
C       Hart, J.F. et al, 'Computer Approximations', Wiley 1968
C
C       Programmer: Alan Miller
C
C       Latest revision - 30 March 1986
C
        IMPLICIT DOUBLE PRECISION (A-H, O-Z)
        DATA P0, P1, P2, P3, P4, P5, P6/220.20 68679 12376 1D0,
     *    221.21 35961 69931 1D0, 112.07 92914 97870 9D0,
     *    33.912 86607 83830 0D0, 6.3739 62203 53165 0D0,
     *    .70038 30644 43688 1D0, .35262 49659 98910 9D-01/,
     *    Q0, Q1, Q2, Q3, Q4, Q5, Q6, Q7/440.41 37358 24752 2D0,
     *    793.82 65125 19948 4D0, 637.33 36333 78831 1D0,
     *    296.56 42487 79673 7D0, 86.780 73220 29460 8D0,
     *    16.064 17757 92069 5D0, 1.7556 67163 18264 2D0,
     *    .88388 34764 83184 4D-1/,
     *    CUTOFF/7.071D0/, ROOT2PI/2.5066 28274 63100 1D0/
C
        ZABS = ABS(Z)
C
C       |Z| > 37.
C
        IF (ZABS .GT. 37.D0) THEN
          PDF = 0.D0
          IF (Z .GT. 0.D0) THEN
            P = 1.D0
            Q = 0.D0
          ELSE
            P = 0.D0
            Q = 1.D0
          END IF
          RETURN
        END IF
C
C       |Z| <= 37.
C
        EXPNTL = EXP(-0.5D0*ZABS**2)
        PDF = EXPNTL/ROOT2PI
C
C       |Z| < CUTOFF = 10/sqrt(2).
C
        IF (ZABS .LT. CUTOFF) THEN
          P = EXPNTL*((((((P6*ZABS + P5)*ZABS + P4)*ZABS + P3)*ZABS +
     *          P2)*ZABS + P1)*ZABS + P0)/(((((((Q7*ZABS + Q6)*ZABS +
     *          Q5)*ZABS + Q4)*ZABS + Q3)*ZABS + Q2)*ZABS + Q1)*ZABS +
     *          Q0)
C
C       |Z| >= CUTOFF.
C
        ELSE
          P = PDF/(ZABS + 1.D0/(ZABS + 2.D0/(ZABS + 3.D0/(ZABS + 4.D0/
     *          (ZABS + 0.65D0)))))
        END IF
C
        IF (Z .LT. 0.D0) THEN
          Q = 1.D0 - P
        ELSE
          Q = P
          P = 1.D0 - Q
        END IF
        RETURN
        END
@}

The only thing to do now is find an efficient way to divide this
up into equal areas.
Symmetry takes care of half of the problem, so we could run along 
our array of wanted points estimating the next value based on
the one before.
Need to just specify the number of points to compute for.

\chapter{Gui interfaces}

For pattern fitting we need 2d plots of x versus y.

For dealing with 2d datasets we need to plot images,
using color or a surface.

For density maps we need to show a 4 d object (x,y,z position
and value).

For crystal structure we need to draw the atoms.

For building models we need to develop ways for the 
user to inspect the model and visualise the restraints
and constraints.

We do not want to tie the linarp code to any particular
gui library.
That means that we need to think carefully about how the 
code is structured.
Our plotting routines, whatever they are, need to just 
pull the information they want to display out of the 
objects we are using.
Similarly, all of our objects should also be available 
from the command line.
The gui could be designed to just display the various objects
and allow you to interact with them, without needing to 
use the python interpreter directly.

If an object contains things which can logically be 
placed on a plot it needs to store them in a way 
that something else can read them.
For peak arrays we have positions and intensities.




\section{2D plots}

There are a wide variety of 2D plotting options for 
python. 
Most developed seems to the the matplotlib package 
from John Hunter.

A two dimensional tk based plotting widget from that package is

@o embedding_in_tk.py
@{
"""
From the matplotlib examples - modified for mouse
"""
import matplotlib
matplotlib.use('TkAgg')

from matplotlib.numerix import arange, sin, cos, pi, searchsorted, sqrt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure

import Tkinter as Tk
import tkFileDialog
import sys,os,time

import powderdata
import profvalplusback
import densemodelfit
import Numeric

class twodplot:
   def __init__(self,dat):
      self.f = Figure(figsize=(5,4), dpi=100)
      self.a = self.f.add_subplot(111)
      self.dat=dat
      self.model=None
      self.a.plot(self.dat.x,self.dat.y)
      if self.dat.d.has_key("xunits"):
         self.a.set_xlabel(dat.d["xunits"])
      if self.dat.d.has_key("yunits"):
         self.a.set_ylabel(dat.d["yunits"])
      if self.dat.d.has_key("title"):
         self.a.set_title(dat.d["title"])
         self.title=dat.d["title"]
      elif self.dat.d.has_key("fromfile"):
         self.a.set_title(dat.d["fromfile"])
         self.title=dat.d["fromfile"]

      # a tk.DrawingArea
      self.canvas = FigureCanvasTkAgg(self.f, master=root)
      self.canvas.show()
      self.tkc=self.canvas.get_tk_widget()
      self.tkc.pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)
      self.tkc.bind("<ButtonPress-1>",self.on_down)
      self.tkc.bind("<ButtonRelease-1>",self.on_up)
      self.tkc.bind("<Button1-Motion>",self.on_move)
      self.tkc.bind("<ButtonPress-2>",self.on_2)
      self.tkc.bind("<ButtonPress-3>",self.on_3)
      self.rubberbandbox=None
      self.label=Tk.Label(master=root, text="Status Bar")
      self.label.pack(side=Tk.TOP)
      self.bf=Tk.Frame(master=root)
      Tk.Button(master=self.bf, text='Quit', command=root.destroy).pack(side=Tk.RIGHT)
      Tk.Button(master=self.bf, text='LogY', command=self.logy).pack(side=Tk.LEFT)
      Tk.Button(master=self.bf, text='LogX', command=self.logx).pack(side=Tk.LEFT)
      Tk.Button(master=self.bf, text='Open xye File', command=self.openxye).pack(side=Tk.LEFT)
      Tk.Button(master=self.bf, text='Open mca File', command=self.openmca).pack(side=Tk.LEFT)
      Tk.Button(master=self.bf, text='Open DILS File', command=self.opendil).pack(side=Tk.LEFT)
      Tk.Button(master=self.bf, text="Fit peak", command=self.estimate).pack(side=Tk.LEFT)
      Tk.Button(master=self.bf, text="Refine prev peak", command=self.pf).pack(side=Tk.LEFT)
      self.bf.pack(side=Tk.BOTTOM)
      self.xfull=self.a.get_xlim()
      self.yfull=self.a.get_ylim()
      self.model=None


   def openxye(self):
      fn=tkFileDialog.askopenfilename(initialdir=os.getcwd())
      self.dat=epffile.epffile(fn)
      self.a.cla()
      time.sleep(1)
      self.a.plot(self.dat.x,self.dat.y)
      self.xfull=self.a.get_xlim()
      self.yfull=self.a.get_ylim()
      self.a.set_xlim(self.xfull)
      self.a.set_ylim(self.yfull)
      self.canvas.show()

   def openmca(self):
      fn=tkFileDialog.askopenfilename(initialdir=os.getcwd())
      self.dat=mcadata.mcadata(fn)
      self.a.cla()
      self.a.plot(self.dat.x,self.dat.y)
      self.xfull=self.a.get_xlim()
      self.yfull=self.a.get_ylim()
      self.a.set_xlim(self.xfull)
      self.a.set_ylim(self.yfull)
      self.canvas.show()

   def opendil(self):
      import solventrefinedata, solventmodel

      ciifile=tkFileDialog.askopenfilename(initialdir=os.getcwd(),title="dils file?")
      cclfile=tkFileDialog.askopenfilename(initialdir=os.getcwd(),title="ccl file?")
      pdbfile=tkFileDialog.askopenfilename(initialdir=os.getcwd(),title="pdb file?")
      self.dat=solventrefinedata.solventrefinedata(ciifile,cclfile,pdbfile)
      # estimate the scale factor:
      calc=sum(self.dat.fsq)
      obs= sum(self.dat.ciiobject.Ihkl)
      s = obs/calc
      self.model = solventmodel.solventmodel(Asolv=0.89,Bsolv=40.0,scale=s)
      self.optimiser=densemodelfit.densemodelfit(model=self.model,data=self.dat)
      self.optimiser.refine(ncycles=10)
      self.a.cla()
      self.a.plot(self.dat.x,self.dat.y)
      self.xfull=self.a.get_xlim()
      self.yfull=self.a.get_ylim()
      self.a.set_xlim(self.xfull)
      self.a.set_ylim(self.yfull)
      self.replot()



   def pf(self):
      self.optimiser.refine(ncycles=1)
      self.replot()

   def estimate(self):
      self.dat.setrange(self.a.get_xlim())
      self.model=profvalplusback.profvalplusback()
      self.optimiser=densemodelfit.densemodelfit(model=self.model,data=self.dat)
      self.optimiser.refine(ncycles=10)
      self.replot()


   def replot(self):
      xr=self.a.get_xlim()
      yr=self.a.get_ylim()
      self.a.clear()
      self.a.set_title(self.title)
      self.a.plot(self.dat.x,self.dat.y)
      if self.model != None:
         self.a.plot(self.dat.x,self.model.ycalc,'g')
      middle=0.5*(yr[0]+yr[1])
         d = self.dat.y-self.model.ycalc
         d = d * self.dat.weightmatvec(d)
         d = d * abs(yr[0]-yr[1]) / Numeric.maximum.reduce(d)
         self.a.plot(self.dat.x, d+middle ,'r')
      self.a.set_xlim(xr)
      self.a.set_ylim(yr)
      self.canvas.show()
      if self.model != None:
         self.label.config(text="Position: %f +/- %f" % (
                   self.model.get_value('center'),
                   self.model.get_errorbar('center'))  )




   def logy(self): 
# FIXME - clip negative values before making logscaled
      if self.a.yaxis.is_log():
         self.a.set_yscale('linear')
      else:
         self.a.set_yscale('log')
      self.canvas.show()

   def logx(self): 
      if self.a.xaxis.is_log():
         self.a.set_xscale('linear')
      else:
         self.a.set_xscale('log')
      self.canvas.show()

   def on_3(self,event):
#      self.dat.setrange()
#      self.a.cla()
#      self.a.plot(self.dat.x, self.dat.y)
      self.a.set_xlim(self.xfull)
      self.a.set_ylim(self.yfull)
#      self.model=None
      self.canvas.show()

   def on_2(self,event):
      height = self.f.bbox.height()
      x, y = event.x, height-event.y
      (xd,yd)= self.a.transData.inverse_xy_tup( (x,y) )
      self.label.config(text="Clicked at x=%f y=%f"%(xd,yd))

   # Callback functions for mouse
   def on_down(self,event):
      # get the x and y coords, flip y from top to bottom
      height = self.f.bbox.height()
      x, y = event.x, height-event.y
      # transData transforms data coords to display coords.  Use the
      # inverse method to transform back
      (self.xd,self.yd)= self.a.transData.inverse_xy_tup( (x,y) )
      # print "print mouse down at", t, val
      # rubber banding:
      if self.rubberbandbox!=None: self.tkc.delete(self.rubberbandbox)
      self.startx=self.tkc.canvasx(event.x)
      self.starty=self.tkc.canvasx(event.y)

   def on_move(self,event):
      x = self.tkc.canvasx(event.x)
      y = self.tkc.canvasy(event.y)
      if (self.startx != event.x)  and (self.starty != event.y) :
         if self.rubberbandbox!=None: 
            self.tkc.delete(self.rubberbandbox)
         self.rubberbandbox = self.tkc.create_rectangle(self.startx, self.starty, x, y, outline='green')
         # this flushes the output, making sure that
         # the rectangle makes it to the screen
         # before the next event is handled
         


   def on_up(self,event):
      # get the x and y coords, flip y from top to bottom
      self.tkc.delete(self.rubberbandbox)
      height = self.f.bbox.height()
      x, y = event.x, height-event.y
      # transData transforms data coords to display coords.  Use the
      # inverse method to transform back
      (self.xu,self.yu) = self.a.transData.inverse_xy_tup( (x,y) )
      if self.xu != self.xd and self.yu != self.yd: 
         # rescale
         xr=[self.xd,self.xu];xr.sort()
         yr=[self.yd,self.yu];yr.sort()
         self.a.set_xlim(xr)
         self.a.set_ylim(yr)
         self.canvas.show()


if __name__=="__main__":
   import epffile, powbase, mcadata, ciidata
   if len(sys.argv)<3:
      print "Usage: %s filename format"%(sys.argv[0])
      x=arange(0.0,3.0,0.01)
      dat=epffile.powderdata(x,sin(2*pi*x)+5,sqrt(sin(2*pi*x)+5),{ "title":"sin x" } )
   else:
      try:
         if sys.argv[2]=="powbase":
            dat=powbase.powbase(sys.argv[1])
         if sys.argv[2]=="epf":
            dat=epffile.epffile(sys.argv[1])
         if sys.argv[2]=="mca":
            dat=mcadata.mcadata(sys.argv[1])
      except:
         print "Could not read your file %s" % (sys.argv[1])
         raise

   root = Tk.Tk()
   root.wm_title("Linarp Is Not A Rietveld Program")
   p=twodplot(dat)

   Tk.mainloop()

@}



\chapter{Technicalities}

Building and maintaining a large project is going to need some
clever stuff for testing. 
At the minimum we need to keep track of the various things to
build, how to build them and how to test them.
Absolutely everything is going to depend on the one master file
in nuweb format.
A makefile will do everything else we want.

@o Makefile -t
@{
all: linarp.dvi linarp.ps linarp.pdf

linarp.pdf : linarp.dvi
	dvipdfm linarp

linarp.ps : linarp.dvi
	dvips linarp

linarp.dvi : linarp.tex
	latex linarp
	nuweb linarp
	latex linarp
	latex linarp

linarp.w : reference.py.out 
	nuweb linarp
	echo "Running python"
	python reference.py > reference.py.out
	echo "Python finished OK"

linarp.tex : linarp.w
	nuweb linarp
@}

The setup script for building extensions.

@o setup.py
@{
import os 
os.system('g77 -c profval.f')              # Nasty hack - needs g77 on path 
os.system('ar -cr libprofval.a profval.o')   # assume libprofval.a link compatible

from distutils.core import setup, Extension

m = Extension("profval", sources = ['cprofval.c'], 
               libraries=['profval'], library_dirs=['.'])


setup(name = 'profval',
      version = '0.1',
      description = "Larry Finger's low angle peakshape",
# Clearly the routine is from Larry Finger and the wrapper is 
# from JPW
      author = 'Jon Wright',
      author_email = 'wright@@esrf.fr',
      ext_modules = [m])

os.system('g77 -c -O2 libchol.f')              # Nasty hack - needs g77 on path 
os.system('ar -cr libchol.a libchol.o')   # assume libprofval.a link compatible


n = Extension("pylibchol", sources = ['pylibchol.c'], 
               libraries=['chol'], library_dirs=['.'])


setup(name = 'pylibchol',
      version = '0.1',
      description = "Sparse Cholesky routines for CII matrices" ,
      author = 'Jon Wright',
      author_email = 'wright@@esrf.fr',
      ext_modules = [n])

os.system("g77 -c bispev.f")
os.system("ar -cr libsplines.a bispev.o")

module1 = Extension("_splines", sources = ['splines.c'], 
                     libraries = ["splines"], library_dirs = ["."] )

setup(name = '_splines',
      version = '0.1',
      description = 'Spline functions',
      ext_modules = [module1])



module2 = Extension("_connectedpixels", sources = ['connectedpixels.c'])

setup(name = 'connectedpixels',
      version = '0.1',
      description = 'connectedpixels',
      ext_modules = [module2])

@}


\section{Fast and accurate timing routine for benchmarking}

Found it on the web, used it when working on the 2D peaksearching.

@O rdtsc.c
@{
/* From: Tom Burgess (Tom_Burgess@@bc.sympatico.ca)
   Subject: DJGPP RDTSC demo (Pentium-only, ~100 lines) 
   Newsgroups: comp.os.msdos.djgpp
   Date: 1997/04/20 

Hi, here's some code that might be useful to some for low-level
Pentium optimization. If you get weird results, look carefully at
what is known to be in cache when the code executes, code & data
alignment, cache line conflicts, AGIs etc. Agner Fog warns that
RDTSC doesn't work with virtual 86 mode but I've noted no problems
with win95 dos shell, RHIDE or whatever. He also points out
special Pentium Pro considerations which I have not addressed.
Check out: http://announce.com/agner/assem/assem.html
 regards, tom  */

#define RDTSC1(dest) \
__asm__(".byte 0x0F, 0x31\n\t"\
        "movl    %%eax, (%%edi)\n\t"\
        "movl    %%edx, 4(%%edi)\n\t"\
        "cld \n\t"\
        "nop \n\t nop \n\t nop \n\t"\
        "nop \n\t nop \n\t nop \n\t"\
 : : "D" (dest) : "eax", "edx")

/* use RDTSC2 immediately after the code under test. The clc 
is a non-pairable filler that also elimate potential shadow effects */

 #define RDTSC2(dest) \
__asm__("clc \n\t"\
    ".byte 0x0F, 0x31\n\t"\
        "movl    %%eax, (%%edi)\n\t"\
        "movl    %%edx, 4(%%edi)\n\t"\
 : : "D" (dest) : "eax", "edx")


/* demo main() - compare sqrt(x) with equivalent pow(x, 0.5) */
/* On P120 get 150 cycles for sqrt, 306 for pow. Strange, since */
/* FSQRT should only be about 70 cycles */

int main()
{
 unsigned long long t1, t2;
 int overhead, ncycles int i;
 double a[10], b[10], c = 123456789.0, d = 999.0, q;
 RDTSC1(&t1); /* just want to get get t1 & t2 into L1 cache */
 RDTSC2(&t2);
 RDTSC1(&t1); /* Measure overhead */
 RDTSC2(&t2);
 overhead = (t2 - t1);
 printf ("\nRDTSC overhead = %d cycles\n", overhead);
 q = sqrt(d); /* get a & sqrt code into cache */
 RDTSC1(&t1); /* Measure sqrt() timing */
 for (i=0;i<10;i++)a[i] = sqrt(c);
 RDTSC2(&t2);
 ncycles = (t2 - t1)/10 - overhead;
 printf ("a = sqrt(c): %d cycles\n", ncycles);
 q = pow(d, 0.5); /* get b & pow code into cache */
 RDTSC1(&t1); /* Measure pow() timing */
 for(i=0;i<10;i++) b[i] = pow(c, 0.5);
 RDTSC2(&t2);
 ncycles = (t2 - t1)/10 - overhead;
 printf ("b = pow(c, 0.5): %d cycles\n", ncycles);
 printf ("a = %f, b = %f\n", a[0], b[0]); /* check results */
  /* and prevent optimizer from eliminating code under test */
 return(0);
}
@}


So for our use we want a header file:
@o rdtsc.h
@{
unsigned long long t1, t2;
int ncycles;
#define RDTSC1(dest) \
__asm__(".byte 0x0F, 0x31\n\t"\
        "movl    %%eax, (%%edi)\n\t"\
        "movl    %%edx, 4(%%edi)\n\t"\
        "cld \n\t"\
        "nop \n\t nop \n\t nop \n\t"\
        "nop \n\t nop \n\t nop \n\t"\
 : : "D" (dest) : "eax", "edx")

#define RDTSC2(dest) \
__asm__("clc \n\t"\
    ".byte 0x0F, 0x31\n\t"\
        "movl    %%eax, (%%edi)\n\t"\
        "movl    %%edx, 4(%%edi)\n\t"\
 : : "D" (dest) : "eax", "edx")

/* usage: RDTSC1(&t1); to start
          RDTSC2(&t2); to stop
          ncycles=t2-t1
*/
@}



\section{License}

Any parts of CCSL which creep into this program were 
copyrighted by the following notice:

\begin{verbatim}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
!                                                                !
!                          DISCLAIMER                            !
!  This software is distributed in the hope that it will be      !
!  useful but \emph{without any warranty}. The author(s) do not  !
!  accept responsibility  to anyone for the consequences of      !
!  using it or for whether it serves  any particular purpose or  !
!  works at all. No warranty is made about  the software or its  !
!  performance.                                                  !
!                                                                !
!                          COPYING                               !
!  Use and copying of this software and the preparation of       !
!  derivative works based on it are permitted, so long as the    !
!  following conditions are met:                                 !
!                                                                !
!    1. The copyright notice and this entire notice are          !
!       included intact and prominently carried on all           !
!       copies and supporting documentation.                     !
!   2.  No fees or compensation are charged for use,             !
!       copies, or access to this software. You may charge       !
!       a nominal distribution fee for the physical act of       !
!       transferring a copy, but you may not charge for the      !
!       programs themselves.                                     !
!    3. If you modify this software, you must cause the          !
!       modified file(s) to carry notices describing the         !
!       changes, who made  the changes, and the date of          !
!       those changes.                                           !
!                                                                !
!                         AUTHORS                                !
!  The CCSL library was written by J.C. Matthewman and           !
!  P.J. Brown with contributions from W.I.F. David, J.B.         !
!  Forsyth and J.H. Matthewman.                                  !
!  Copyright c  1999--2000. All rights reserved.                 !
!                                                                !             
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
\end{verbatim}

I think that the main thing which is likely to come into 
linarp from CCSL is contained the in the subroutine matdil.f (and
it's friends) which was written entirely from scratch by the 
author.

It is the authors firm belief that scientific software is utterly
useless unless you can find out \emph{exactly} how it works.
As such, you actually need to be able to read the source code
of a program which does scientific data analysis, if you want
to know what the program does.
Manuals only tell you what the author wanted the program to do,
you need to know what it does and how it does it if you 
want to interpret any results.
For better or worse, the GPL license was chosen to ensure that the
source code should follow the program around. That seems to mean
that any file created needs a copyright notice attached to
it, so versions of that notice for c, fortran and python are
defined here.


@d fortrancopyright
@{
!  Copyright (C) 2004 by Jonathan Wright
!     Grenoble, France
!     email: wright@@esrf.fr
!
!   This file is part of Linarp.
!
!   Linarp is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation; either version 2 of the License, or
!   (at your option) any later version.
!
!   Linarp is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with Linarp; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
@}

@d ccopyright
@{
/* 
@< fortrancopyright @>
*/
@}

@d pycopyright
@{
"""
@< fortrancopyright @>
"""
@}

Now for the GPL itself. Rather long, isn't it.

@O LICENSE
@{
		    GNU GENERAL PUBLIC LICENSE
		       Version 2, June 1991

 Copyright (C) 1989, 1991 Free Software Foundation, Inc.
                       59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 Everyone is permitted to copy and distribute verbatim copies
 of this license document, but changing it is not allowed.

			    Preamble

  The licenses for most software are designed to take away your
freedom to share and change it.  By contrast, the GNU General Public
License is intended to guarantee your freedom to share and change free
software--to make sure the software is free for all its users.  This
General Public License applies to most of the Free Software
Foundation's software and to any other program whose authors commit to
using it.  (Some other Free Software Foundation software is covered by
the GNU Library General Public License instead.)  You can apply it to
your programs, too.

  When we speak of free software, we are referring to freedom, not
price.  Our General Public Licenses are designed to make sure that you
have the freedom to distribute copies of free software (and charge for
this service if you wish), that you receive source code or can get it
if you want it, that you can change the software or use pieces of it
in new free programs; and that you know you can do these things.

  To protect your rights, we need to make restrictions that forbid
anyone to deny you these rights or to ask you to surrender the rights.
These restrictions translate to certain responsibilities for you if you
distribute copies of the software, or if you modify it.

  For example, if you distribute copies of such a program, whether
gratis or for a fee, you must give the recipients all the rights that
you have.  You must make sure that they, too, receive or can get the
source code.  And you must show them these terms so they know their
rights.

  We protect your rights with two steps: (1) copyright the software, and
(2) offer you this license which gives you legal permission to copy,
distribute and/or modify the software.

  Also, for each author's protection and ours, we want to make certain
that everyone understands that there is no warranty for this free
software.  If the software is modified by someone else and passed on, we
want its recipients to know that what they have is not the original, so
that any problems introduced by others will not reflect on the original
authors' reputations.

  Finally, any free program is threatened constantly by software
patents.  We wish to avoid the danger that redistributors of a free
program will individually obtain patent licenses, in effect making the
program proprietary.  To prevent this, we have made it clear that any
patent must be licensed for everyone's free use or not licensed at all.

  The precise terms and conditions for copying, distribution and
modification follow.

		    GNU GENERAL PUBLIC LICENSE
   TERMS AND CONDITIONS FOR COPYING, DISTRIBUTION AND MODIFICATION

  0. This License applies to any program or other work which contains
a notice placed by the copyright holder saying it may be distributed
under the terms of this General Public License.  The "Program", below,
refers to any such program or work, and a "work based on the Program"
means either the Program or any derivative work under copyright law:
that is to say, a work containing the Program or a portion of it,
either verbatim or with modifications and/or translated into another
language.  (Hereinafter, translation is included without limitation in
the term "modification".)  Each licensee is addressed as "you".

Activities other than copying, distribution and modification are not
covered by this License; they are outside its scope.  The act of
running the Program is not restricted, and the output from the Program
is covered only if its contents constitute a work based on the
Program (independent of having been made by running the Program).
Whether that is true depends on what the Program does.

  1. You may copy and distribute verbatim copies of the Program's
source code as you receive it, in any medium, provided that you
conspicuously and appropriately publish on each copy an appropriate
copyright notice and disclaimer of warranty; keep intact all the
notices that refer to this License and to the absence of any warranty;
and give any other recipients of the Program a copy of this License
along with the Program.

You may charge a fee for the physical act of transferring a copy, and
you may at your option offer warranty protection in exchange for a fee.

  2. You may modify your copy or copies of the Program or any portion
of it, thus forming a work based on the Program, and copy and
distribute such modifications or work under the terms of Section 1
above, provided that you also meet all of these conditions:

    a) You must cause the modified files to carry prominent notices
    stating that you changed the files and the date of any change.

    b) You must cause any work that you distribute or publish, that in
    whole or in part contains or is derived from the Program or any
    part thereof, to be licensed as a whole at no charge to all third
    parties under the terms of this License.

    c) If the modified program normally reads commands interactively
    when run, you must cause it, when started running for such
    interactive use in the most ordinary way, to print or display an
    announcement including an appropriate copyright notice and a
    notice that there is no warranty (or else, saying that you provide
    a warranty) and that users may redistribute the program under
    these conditions, and telling the user how to view a copy of this
    License.  (Exception: if the Program itself is interactive but
    does not normally print such an announcement, your work based on
    the Program is not required to print an announcement.)

These requirements apply to the modified work as a whole.  If
identifiable sections of that work are not derived from the Program,
and can be reasonably considered independent and separate works in
themselves, then this License, and its terms, do not apply to those
sections when you distribute them as separate works.  But when you
distribute the same sections as part of a whole which is a work based
on the Program, the distribution of the whole must be on the terms of
this License, whose permissions for other licensees extend to the
entire whole, and thus to each and every part regardless of who wrote it.

Thus, it is not the intent of this section to claim rights or contest
your rights to work written entirely by you; rather, the intent is to
exercise the right to control the distribution of derivative or
collective works based on the Program.

In addition, mere aggregation of another work not based on the Program
with the Program (or with a work based on the Program) on a volume of
a storage or distribution medium does not bring the other work under
the scope of this License.

  3. You may copy and distribute the Program (or a work based on it,
under Section 2) in object code or executable form under the terms of
Sections 1 and 2 above provided that you also do one of the following:

    a) Accompany it with the complete corresponding machine-readable
    source code, which must be distributed under the terms of Sections
    1 and 2 above on a medium customarily used for software interchange; or,

    b) Accompany it with a written offer, valid for at least three
    years, to give any third party, for a charge no more than your
    cost of physically performing source distribution, a complete
    machine-readable copy of the corresponding source code, to be
    distributed under the terms of Sections 1 and 2 above on a medium
    customarily used for software interchange; or,

    c) Accompany it with the information you received as to the offer
    to distribute corresponding source code.  (This alternative is
    allowed only for noncommercial distribution and only if you
    received the program in object code or executable form with such
    an offer, in accord with Subsection b above.)

The source code for a work means the preferred form of the work for
making modifications to it.  For an executable work, complete source
code means all the source code for all modules it contains, plus any
associated interface definition files, plus the scripts used to
control compilation and installation of the executable.  However, as a
special exception, the source code distributed need not include
anything that is normally distributed (in either source or binary
form) with the major components (compiler, kernel, and so on) of the
operating system on which the executable runs, unless that component
itself accompanies the executable.

If distribution of executable or object code is made by offering
access to copy from a designated place, then offering equivalent
access to copy the source code from the same place counts as
distribution of the source code, even though third parties are not
compelled to copy the source along with the object code.

  4. You may not copy, modify, sublicense, or distribute the Program
except as expressly provided under this License.  Any attempt
otherwise to copy, modify, sublicense or distribute the Program is
void, and will automatically terminate your rights under this License.
However, parties who have received copies, or rights, from you under
this License will not have their licenses terminated so long as such
parties remain in full compliance.

  5. You are not required to accept this License, since you have not
signed it.  However, nothing else grants you permission to modify or
distribute the Program or its derivative works.  These actions are
prohibited by law if you do not accept this License.  Therefore, by
modifying or distributing the Program (or any work based on the
Program), you indicate your acceptance of this License to do so, and
all its terms and conditions for copying, distributing or modifying
the Program or works based on it.

  6. Each time you redistribute the Program (or any work based on the
Program), the recipient automatically receives a license from the
original licensor to copy, distribute or modify the Program subject to
these terms and conditions.  You may not impose any further
restrictions on the recipients' exercise of the rights granted herein.
You are not responsible for enforcing compliance by third parties to
this License.

  7. If, as a consequence of a court judgment or allegation of patent
infringement or for any other reason (not limited to patent issues),
conditions are imposed on you (whether by court order, agreement or
otherwise) that contradict the conditions of this License, they do not
excuse you from the conditions of this License.  If you cannot
distribute so as to satisfy simultaneously your obligations under this
License and any other pertinent obligations, then as a consequence you
may not distribute the Program at all.  For example, if a patent
license would not permit royalty-free redistribution of the Program by
all those who receive copies directly or indirectly through you, then
the only way you could satisfy both it and this License would be to
refrain entirely from distribution of the Program.

If any portion of this section is held invalid or unenforceable under
any particular circumstance, the balance of the section is intended to
apply and the section as a whole is intended to apply in other
circumstances.

It is not the purpose of this section to induce you to infringe any
patents or other property right claims or to contest validity of any
such claims; this section has the sole purpose of protecting the
integrity of the free software distribution system, which is
implemented by public license practices.  Many people have made
generous contributions to the wide range of software distributed
through that system in reliance on consistent application of that
system; it is up to the author/donor to decide if he or she is willing
to distribute software through any other system and a licensee cannot
impose that choice.

This section is intended to make thoroughly clear what is believed to
be a consequence of the rest of this License.

  8. If the distribution and/or use of the Program is restricted in
certain countries either by patents or by copyrighted interfaces, the
original copyright holder who places the Program under this License
may add an explicit geographical distribution limitation excluding
those countries, so that distribution is permitted only in or among
countries not thus excluded.  In such case, this License incorporates
the limitation as if written in the body of this License.

  9. The Free Software Foundation may publish revised and/or new versions
of the General Public License from time to time.  Such new versions will
be similar in spirit to the present version, but may differ in detail to
address new problems or concerns.

Each version is given a distinguishing version number.  If the Program
specifies a version number of this License which applies to it and "any
later version", you have the option of following the terms and conditions
either of that version or of any later version published by the Free
Software Foundation.  If the Program does not specify a version number of
this License, you may choose any version ever published by the Free Software
Foundation.

  10. If you wish to incorporate parts of the Program into other free
programs whose distribution conditions are different, write to the author
to ask for permission.  For software which is copyrighted by the Free
Software Foundation, write to the Free Software Foundation; we sometimes
make exceptions for this.  Our decision will be guided by the two goals
of preserving the free status of all derivatives of our free software and
of promoting the sharing and reuse of software generally.

			    NO WARRANTY

  11. BECAUSE THE PROGRAM IS LICENSED FREE OF CHARGE, THERE IS NO WARRANTY
FOR THE PROGRAM, TO THE EXTENT PERMITTED BY APPLICABLE LAW.  EXCEPT WHEN
OTHERWISE STATED IN WRITING THE COPYRIGHT HOLDERS AND/OR OTHER PARTIES
PROVIDE THE PROGRAM "AS IS" WITHOUT WARRANTY OF ANY KIND, EITHER EXPRESSED
OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.  THE ENTIRE RISK AS
TO THE QUALITY AND PERFORMANCE OF THE PROGRAM IS WITH YOU.  SHOULD THE
PROGRAM PROVE DEFECTIVE, YOU ASSUME THE COST OF ALL NECESSARY SERVICING,
REPAIR OR CORRECTION.

  12. IN NO EVENT UNLESS REQUIRED BY APPLICABLE LAW OR AGREED TO IN WRITING
WILL ANY COPYRIGHT HOLDER, OR ANY OTHER PARTY WHO MAY MODIFY AND/OR
REDISTRIBUTE THE PROGRAM AS PERMITTED ABOVE, BE LIABLE TO YOU FOR DAMAGES,
INCLUDING ANY GENERAL, SPECIAL, INCIDENTAL OR CONSEQUENTIAL DAMAGES ARISING
OUT OF THE USE OR INABILITY TO USE THE PROGRAM (INCLUDING BUT NOT LIMITED
TO LOSS OF DATA OR DATA BEING RENDERED INACCURATE OR LOSSES SUSTAINED BY
YOU OR THIRD PARTIES OR A FAILURE OF THE PROGRAM TO OPERATE WITH ANY OTHER
PROGRAMS), EVEN IF SUCH HOLDER OR OTHER PARTY HAS BEEN ADVISED OF THE
POSSIBILITY OF SUCH DAMAGES.

		     END OF TERMS AND CONDITIONS

	    How to Apply These Terms to Your New Programs

  If you develop a new program, and you want it to be of the greatest
possible use to the public, the best way to achieve this is to make it
free software which everyone can redistribute and change under these terms.

  To do so, attach the following notices to the program.  It is safest
to attach them to the start of each source file to most effectively
convey the exclusion of warranty; and each file should have at least
the "copyright" line and a pointer to where the full notice is found.

    <one line to give the program's name and a brief idea of what it does.>
    Copyright (C) <year>  <name of author>

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


Also add information on how to contact you by electronic and paper mail.

If the program is interactive, make it output a short notice like this
when it starts in an interactive mode:

    Gnomovision version 69, Copyright (C) year name of author
    Gnomovision comes with ABSOLUTELY NO WARRANTY; for details type `show w'.
    This is free software, and you are welcome to redistribute it
    under certain conditions; type `show c' for details.

The hypothetical commands `show w' and `show c' should show the appropriate
parts of the General Public License.  Of course, the commands you use may
be called something other than `show w' and `show c'; they could even be
mouse-clicks or menu items--whatever suits your program.

You should also get your employer (if you work as a programmer) or your
school, if any, to sign a "copyright disclaimer" for the program, if
necessary.  Here is a sample; alter the names:

  Yoyodyne, Inc., hereby disclaims all copyright interest in the program
  `Gnomovision' (which makes passes at compilers) written by James Hacker.

  <signature of Ty Coon>, 1 April 1989
  Ty Coon, President of Vice

This General Public License does not permit incorporating your program into
proprietary programs.  If your program is a subroutine library, you may
consider it more useful to permit linking proprietary applications with the
library.  If this is what you want to do, use the GNU Library General
Public License instead of this License.
@}

\section{Acknowledgements}

This work was carried out while the author was employed by the 
European Synchrotron Radiation Facility. 
The sparse matrix decompositions were written while working
and BM16 and a great deal of useful advice, encouragement and
support from Andy Fitch is gratefully acknowledged.
The routines for actually carrying out the sparse cholesky 
decomposition are remarkably similar to those in the 
book by George and Liu~\cite{GeorgeandLiu}, and they were 
written after reading the material in the early chapters. 
The actual implementation is slightly different and a more 
awkward array indexing was used as unfortunately I did
not notice the code until it was too late!

Much of the methodology has been inspired by the work of the 
ISIS/RAL group (David, Sivia and Shankland) 
on the use of Pawley refinements.
Bill David originally suggested to the author that a sparse
Pawley implementation would be worthwhile, and also showed
me an unpublished manuscript demonstrating the relationship
between the correlated integrated intensity $\chi^2$ and
the Rietveld $\chi^2$.

\appendix{Profval tables}

This should intentionally overrun a page when printing.
Replace @@d with @@D if you really want all the pages
of numbers.

@d profvaldata
@{
!
! Values for the abscissas and weights of the Gauss-Legendre
!  N-point quadrature formula have been precomputed using routine
!  Gauleg from "Numerical Recipes" (Press, Flannery, Teukolsky
!  and Vetterling, 1986, Cambridge University Press,
!  ISBN 0 521 30811 9), and are stored in the DATA statements
!  for XP and WP below.
!
      data (xp(i),i=   1,  40)/
     $.2386192E+00,.6612094E+00,.9324695E+00,.1488743E+00,.4333954E+00,
     $.6794096E+00,.8650634E+00,.9739065E+00,.7652652E-01,.2277859E+00,
     $.3737061E+00,.5108670E+00,.6360537E+00,.7463319E+00,.8391170E+00,
     $.9122344E+00,.9639719E+00,.9931286E+00,.3877242E-01,.1160841E+00,
     $.1926976E+00,.2681522E+00,.3419941E+00,.4137792E+00,.4830758E+00,
     $.5494671E+00,.6125539E+00,.6719567E+00,.7273183E+00,.7783057E+00,
     $.8246122E+00,.8659595E+00,.9020988E+00,.9328128E+00,.9579168E+00,
     $.9772599E+00,.9907262E+00,.9982377E+00,.2595977E-01,.7780933E-01/
      data (xp(i),i=  41,  80)/
     $.1294491E+00,.1807400E+00,.2315436E+00,.2817229E+00,.3311428E+00,
     $.3796701E+00,.4271737E+00,.4735258E+00,.5186014E+00,.5622789E+00,
     $.6044406E+00,.6449728E+00,.6837663E+00,.7207165E+00,.7557238E+00,
     $.7886937E+00,.8195375E+00,.8481720E+00,.8745199E+00,.8985103E+00,
     $.9200785E+00,.9391663E+00,.9557223E+00,.9697018E+00,.9810672E+00,
     $.9897879E+00,.9958405E+00,.9992101E+00,.1951138E-01,.5850444E-01,
     $.9740840E-01,.1361640E+00,.1747123E+00,.2129945E+00,.2509524E+00,
     $.2885281E+00,.3256644E+00,.3623048E+00,.3983934E+00,.4338754E+00/
      data (xp(i),i=  81, 120)/
     $.4686966E+00,.5028041E+00,.5361459E+00,.5686713E+00,.6003306E+00,
     $.6310758E+00,.6608599E+00,.6896376E+00,.7173652E+00,.7440003E+00,
     $.7695024E+00,.7938327E+00,.8169541E+00,.8388315E+00,.8594314E+00,
     $.8787226E+00,.8966756E+00,.9132631E+00,.9284599E+00,.9422428E+00,
     $.9545908E+00,.9654851E+00,.9749091E+00,.9828486E+00,.9892913E+00,
     $.9942275E+00,.9976499E+00,.9995538E+00,.1562898E-01,.4687168E-01,
     $.7806858E-01,.1091892E+00,.1402031E+00,.1710801E+00,.2017899E+00,
     $.2323025E+00,.2625881E+00,.2926172E+00,.3223603E+00,.3517885E+00/
      data (xp(i),i= 121, 160)/
     $.3808730E+00,.4095853E+00,.4378974E+00,.4657816E+00,.4932108E+00,
     $.5201580E+00,.5465970E+00,.5725019E+00,.5978475E+00,.6226089E+00,
     $.6467619E+00,.6702830E+00,.6931492E+00,.7153381E+00,.7368281E+00,
     $.7575981E+00,.7776279E+00,.7968979E+00,.8153892E+00,.8330839E+00,
     $.8499645E+00,.8660147E+00,.8812187E+00,.8955616E+00,.9090296E+00,
     $.9216093E+00,.9332885E+00,.9440559E+00,.9539008E+00,.9628137E+00,
     $.9707858E+00,.9778094E+00,.9838775E+00,.9889844E+00,.9931249E+00,
     $.9962951E+00,.9984920E+00,.9997137E+00,.1043694E-01,.3130627E-01/
      data (xp(i),i= 161, 200)/
     $.5216195E-01,.7299491E-01,.9379607E-01,.1145563E+00,.1352667E+00,
     $.1559181E+00,.1765016E+00,.1970082E+00,.2174290E+00,.2377550E+00,
     $.2579774E+00,.2780874E+00,.2980762E+00,.3179352E+00,.3376556E+00,
     $.3572289E+00,.3766466E+00,.3959001E+00,.4149811E+00,.4338813E+00,
     $.4525925E+00,.4711065E+00,.4894151E+00,.5075106E+00,.5253849E+00,
     $.5430303E+00,.5604390E+00,.5776036E+00,.5945165E+00,.6111703E+00,
     $.6275579E+00,.6436720E+00,.6595056E+00,.6750519E+00,.6903041E+00,
     $.7052554E+00,.7198995E+00,.7342299E+00,.7482404E+00,.7619248E+00/
      data (xp(i),i= 201, 240)/
     $.7752773E+00,.7882919E+00,.8009631E+00,.8132853E+00,.8252531E+00,
     $.8368613E+00,.8481049E+00,.8589789E+00,.8694787E+00,.8795996E+00,
     $.8893372E+00,.8986874E+00,.9076460E+00,.9162090E+00,.9243729E+00,
     $.9321340E+00,.9394890E+00,.9464346E+00,.9529678E+00,.9590857E+00,
     $.9647858E+00,.9700655E+00,.9749225E+00,.9793548E+00,.9833603E+00,
     $.9869373E+00,.9900843E+00,.9927999E+00,.9950829E+00,.9969323E+00,
     $.9983473E+00,.9993274E+00,.9998723E+00,.7834291E-02,.2350095E-01,
     $.3916184E-01,.5481311E-01,.7045093E-01,.8607145E-01,.1016708E+00/
      data (xp(i),i= 241, 280)/
     $.1172453E+00,.1327909E+00,.1483040E+00,.1637806E+00,.1792170E+00,
     $.1946095E+00,.2099541E+00,.2252472E+00,.2404850E+00,.2556638E+00,
     $.2707798E+00,.2858293E+00,.3008086E+00,.3157141E+00,.3305421E+00,
     $.3452890E+00,.3599510E+00,.3745247E+00,.3890065E+00,.4033927E+00,
     $.4176799E+00,.4318646E+00,.4459432E+00,.4599124E+00,.4737686E+00,
     $.4875086E+00,.5011288E+00,.5146260E+00,.5279969E+00,.5412382E+00,
     $.5543465E+00,.5673188E+00,.5801518E+00,.5928424E+00,.6053874E+00,
     $.6177838E+00,.6300285E+00,.6421185E+00,.6540509E+00,.6658228E+00/
      data (xp(i),i= 281, 320)/
     $.6774311E+00,.6888732E+00,.7001461E+00,.7112472E+00,.7221736E+00,
     $.7329227E+00,.7434919E+00,.7538786E+00,.7640801E+00,.7740941E+00,
     $.7839181E+00,.7935496E+00,.8029862E+00,.8122257E+00,.8212659E+00,
     $.8301044E+00,.8387391E+00,.8471679E+00,.8553887E+00,.8633995E+00,
     $.8711983E+00,.8787832E+00,.8861524E+00,.8933041E+00,.9002364E+00,
     $.9069477E+00,.9134364E+00,.9197008E+00,.9257394E+00,.9315507E+00,
     $.9371333E+00,.9424859E+00,.9476071E+00,.9524956E+00,.9571503E+00,
     $.9615700E+00,.9657536E+00,.9697002E+00,.9734086E+00,.9768781E+00/
      data (xp(i),i= 321, 360)/
     $.9801078E+00,.9830968E+00,.9858445E+00,.9883502E+00,.9906132E+00,
     $.9926330E+00,.9944091E+00,.9959410E+00,.9972285E+00,.9982712E+00,
     $.9990687E+00,.9996210E+00,.9999281E+00,.5227245E-02,.1568116E-01,
     $.2613337E-01,.3658271E-01,.4702806E-01,.5746827E-01,.6790220E-01,
     $.7832871E-01,.8874665E-01,.9915490E-01,.1095523E+00,.1199377E+00,
     $.1303101E+00,.1406682E+00,.1510109E+00,.1613371E+00,.1716456E+00,
     $.1819354E+00,.1922054E+00,.2024543E+00,.2126811E+00,.2228846E+00,
     $.2330638E+00,.2432175E+00,.2533446E+00,.2634441E+00,.2735147E+00/
      data (xp(i),i= 361, 400)/
     $.2835555E+00,.2935652E+00,.3035429E+00,.3134874E+00,.3233976E+00,
     $.3332725E+00,.3431110E+00,.3529120E+00,.3626744E+00,.3723971E+00,
     $.3820792E+00,.3917194E+00,.4013169E+00,.4108705E+00,.4203792E+00,
     $.4298420E+00,.4392578E+00,.4486255E+00,.4579443E+00,.4672130E+00,
     $.4764306E+00,.4855961E+00,.4947086E+00,.5037670E+00,.5127704E+00,
     $.5217177E+00,.5306079E+00,.5394402E+00,.5482135E+00,.5569269E+00,
     $.5655795E+00,.5741702E+00,.5826982E+00,.5911624E+00,.5995621E+00,
     $.6078963E+00,.6161639E+00,.6243643E+00,.6324964E+00,.6405594E+00/
      data (xp(i),i= 401, 440)/
     $.6485524E+00,.6564744E+00,.6643248E+00,.6721025E+00,.6798068E+00,
     $.6874367E+00,.6949916E+00,.7024704E+00,.7098725E+00,.7171970E+00,
     $.7244432E+00,.7316101E+00,.7386971E+00,.7457033E+00,.7526281E+00,
     $.7594705E+00,.7662300E+00,.7729057E+00,.7794970E+00,.7860030E+00,
     $.7924232E+00,.7987567E+00,.8050030E+00,.8111612E+00,.8172308E+00,
     $.8232111E+00,.8291014E+00,.8349011E+00,.8406095E+00,.8462260E+00,
     $.8517501E+00,.8571811E+00,.8625184E+00,.8677614E+00,.8729095E+00,
     $.8779623E+00,.8829191E+00,.8877794E+00,.8925427E+00,.8972084E+00/
      data (xp(i),i= 441, 480)/
     $.9017761E+00,.9062452E+00,.9106152E+00,.9148857E+00,.9190563E+00,
     $.9231263E+00,.9270955E+00,.9309634E+00,.9347295E+00,.9383934E+00,
     $.9419548E+00,.9454132E+00,.9487683E+00,.9520197E+00,.9551671E+00,
     $.9582100E+00,.9611482E+00,.9639814E+00,.9667092E+00,.9693313E+00,
     $.9718476E+00,.9742575E+00,.9765610E+00,.9787578E+00,.9808476E+00,
     $.9828302E+00,.9847054E+00,.9864729E+00,.9881326E+00,.9896844E+00,
     $.9911279E+00,.9924632E+00,.9936899E+00,.9948081E+00,.9958175E+00,
     $.9967181E+00,.9975097E+00,.9981923E+00,.9987659E+00,.9992302E+00/
      data (xp(i),i= 481, 520)/
     $.9995854E+00,.9998313E+00,.9999680E+00,.3922075E-02,.1176598E-01,
     $.1960917E-01,.2745115E-01,.3529144E-01,.4312955E-01,.5096502E-01,
     $.5879735E-01,.6662606E-01,.7445067E-01,.8227070E-01,.9008566E-01,
     $.9789509E-01,.1056985E+00,.1134954E+00,.1212853E+00,.1290678E+00,
     $.1368423E+00,.1446083E+00,.1523655E+00,.1601134E+00,.1678513E+00,
     $.1755790E+00,.1832958E+00,.1910013E+00,.1986951E+00,.2063767E+00,
     $.2140456E+00,.2217013E+00,.2293434E+00,.2369713E+00,.2445847E+00,
     $.2521830E+00,.2597658E+00,.2673327E+00,.2748830E+00,.2824165E+00/
      data (xp(i),i= 521, 560)/
     $.2899326E+00,.2974308E+00,.3049108E+00,.3123719E+00,.3198139E+00,
     $.3272362E+00,.3346383E+00,.3420199E+00,.3493804E+00,.3567194E+00,
     $.3640365E+00,.3713311E+00,.3786029E+00,.3858515E+00,.3930762E+00,
     $.4002768E+00,.4074528E+00,.4146037E+00,.4217291E+00,.4288285E+00,
     $.4359016E+00,.4429478E+00,.4499667E+00,.4569580E+00,.4639212E+00,
     $.4708558E+00,.4777615E+00,.4846377E+00,.4914841E+00,.4983003E+00,
     $.5050859E+00,.5118403E+00,.5185633E+00,.5252543E+00,.5319131E+00,
     $.5385391E+00,.5451319E+00,.5516912E+00,.5582166E+00,.5647076E+00/
      data (xp(i),i= 561, 600)/
     $.5711639E+00,.5775851E+00,.5839707E+00,.5903203E+00,.5966337E+00,
     $.6029103E+00,.6091498E+00,.6153519E+00,.6215161E+00,.6276420E+00,
     $.6337293E+00,.6397777E+00,.6457866E+00,.6517559E+00,.6576850E+00,
     $.6635737E+00,.6694215E+00,.6752281E+00,.6809932E+00,.6867164E+00,
     $.6923974E+00,.6980357E+00,.7036311E+00,.7091832E+00,.7146916E+00,
     $.7201561E+00,.7255763E+00,.7309518E+00,.7362823E+00,.7415676E+00,
     $.7468072E+00,.7520008E+00,.7571482E+00,.7622490E+00,.7673029E+00,
     $.7723096E+00,.7772688E+00,.7821801E+00,.7870433E+00,.7918581E+00/
      data (xp(i),i= 601, 640)/
     $.7966241E+00,.8013412E+00,.8060089E+00,.8106271E+00,.8151953E+00,
     $.8197134E+00,.8241811E+00,.8285980E+00,.8329640E+00,.8372787E+00,
     $.8415419E+00,.8457533E+00,.8499127E+00,.8540198E+00,.8580743E+00,
     $.8620760E+00,.8660247E+00,.8699201E+00,.8737620E+00,.8775501E+00,
     $.8812842E+00,.8849641E+00,.8885896E+00,.8921603E+00,.8956762E+00,
     $.8991369E+00,.9025424E+00,.9058923E+00,.9091864E+00,.9124246E+00,
     $.9156067E+00,.9187324E+00,.9218016E+00,.9248141E+00,.9277697E+00,
     $.9306682E+00,.9335094E+00,.9362932E+00,.9390194E+00,.9416878E+00/
      data (xp(i),i= 641, 680)/
     $.9442982E+00,.9468506E+00,.9493447E+00,.9517803E+00,.9541574E+00,
     $.9564759E+00,.9587354E+00,.9609360E+00,.9630774E+00,.9651596E+00,
     $.9671823E+00,.9691456E+00,.9710493E+00,.9728932E+00,.9746772E+00,
     $.9764012E+00,.9780652E+00,.9796690E+00,.9812125E+00,.9826957E+00,
     $.9841183E+00,.9854805E+00,.9867820E+00,.9880227E+00,.9892027E+00,
     $.9903218E+00,.9913800E+00,.9923771E+00,.9933133E+00,.9941882E+00,
     $.9950021E+00,.9957547E+00,.9964460E+00,.9970760E+00,.9976447E+00,
     $.9981519E+00,.9985978E+00,.9989822E+00,.9993052E+00,.9995666E+00/
      data (xp(i),i= 681, 720)/
     $.9997666E+00,.9999050E+00,.9999820E+00,.2615810E-02,.7847359E-02,
     $.1307869E-01,.1830967E-01,.2354014E-01,.2876997E-01,.3399902E-01,
     $.3922713E-01,.4445417E-01,.4967999E-01,.5490445E-01,.6012741E-01,
     $.6534873E-01,.7056825E-01,.7578585E-01,.8100137E-01,.8621467E-01,
     $.9142561E-01,.9663405E-01,.1018398E+00,.1070429E+00,.1122429E+00,
     $.1174399E+00,.1226337E+00,.1278242E+00,.1330111E+00,.1381944E+00,
     $.1433739E+00,.1485495E+00,.1537210E+00,.1588884E+00,.1640513E+00,
     $.1692098E+00,.1743636E+00,.1795127E+00,.1846569E+00,.1897960E+00/
      data (xp(i),i= 721, 760)/
     $.1949299E+00,.2000585E+00,.2051816E+00,.2102991E+00,.2154108E+00,
     $.2205166E+00,.2256164E+00,.2307101E+00,.2357974E+00,.2408782E+00,
     $.2459525E+00,.2510200E+00,.2560807E+00,.2611343E+00,.2661808E+00,
     $.2712201E+00,.2762519E+00,.2812761E+00,.2862926E+00,.2913013E+00,
     $.2963021E+00,.3012947E+00,.3062790E+00,.3112550E+00,.3162225E+00,
     $.3211813E+00,.3261313E+00,.3310724E+00,.3360045E+00,.3409273E+00,
     $.3458408E+00,.3507449E+00,.3556393E+00,.3605240E+00,.3653989E+00,
     $.3702637E+00,.3751184E+00,.3799629E+00,.3847969E+00,.3896204E+00/
      data (xp(i),i= 761, 800)/
     $.3944333E+00,.3992353E+00,.4040264E+00,.4088065E+00,.4135754E+00,
     $.4183329E+00,.4230790E+00,.4278136E+00,.4325364E+00,.4372474E+00,
     $.4419464E+00,.4466333E+00,.4513080E+00,.4559703E+00,.4606202E+00,
     $.4652574E+00,.4698819E+00,.4744936E+00,.4790923E+00,.4836778E+00,
     $.4882502E+00,.4928091E+00,.4973546E+00,.5018864E+00,.5064046E+00,
     $.5109088E+00,.5153991E+00,.5198753E+00,.5243372E+00,.5287848E+00,
     $.5332179E+00,.5376364E+00,.5420402E+00,.5464292E+00,.5508032E+00,
     $.5551622E+00,.5595059E+00,.5638343E+00,.5681473E+00,.5724448E+00/
      data (xp(i),i= 801, 840)/
     $.5767266E+00,.5809926E+00,.5852427E+00,.5894768E+00,.5936947E+00,
     $.5978964E+00,.6020817E+00,.6062506E+00,.6104028E+00,.6145384E+00,
     $.6186571E+00,.6227589E+00,.6268437E+00,.6309113E+00,.6349616E+00,
     $.6389945E+00,.6430100E+00,.6470079E+00,.6509880E+00,.6549504E+00,
     $.6588948E+00,.6628211E+00,.6667294E+00,.6706194E+00,.6744910E+00,
     $.6783442E+00,.6821788E+00,.6859947E+00,.6897919E+00,.6935702E+00,
     $.6973295E+00,.7010697E+00,.7047907E+00,.7084924E+00,.7121748E+00,
     $.7158376E+00,.7194809E+00,.7231044E+00,.7267082E+00,.7302921E+00/
      data (xp(i),i= 841, 880)/
     $.7338560E+00,.7373998E+00,.7409234E+00,.7444268E+00,.7479097E+00,
     $.7513722E+00,.7548142E+00,.7582355E+00,.7616360E+00,.7650157E+00,
     $.7683744E+00,.7717121E+00,.7750287E+00,.7783241E+00,.7815982E+00,
     $.7848508E+00,.7880821E+00,.7912917E+00,.7944797E+00,.7976459E+00,
     $.8007903E+00,.8039128E+00,.8070132E+00,.8100916E+00,.8131479E+00,
     $.8161818E+00,.8191934E+00,.8221826E+00,.8251493E+00,.8280935E+00,
     $.8310149E+00,.8339136E+00,.8367895E+00,.8396425E+00,.8424725E+00,
     $.8452794E+00,.8480632E+00,.8508238E+00,.8535611E+00,.8562750E+00/
      data (xp(i),i= 881, 920)/
     $.8589656E+00,.8616325E+00,.8642760E+00,.8668957E+00,.8694918E+00,
     $.8720640E+00,.8746124E+00,.8771368E+00,.8796372E+00,.8821136E+00,
     $.8845658E+00,.8869937E+00,.8893975E+00,.8917768E+00,.8941318E+00,
     $.8964623E+00,.8987683E+00,.9010496E+00,.9033063E+00,.9055383E+00,
     $.9077455E+00,.9099278E+00,.9120852E+00,.9142177E+00,.9163252E+00,
     $.9184075E+00,.9204648E+00,.9224968E+00,.9245036E+00,.9264851E+00,
     $.9284412E+00,.9303720E+00,.9322772E+00,.9341570E+00,.9360111E+00,
     $.9378397E+00,.9396426E+00,.9414198E+00,.9431712E+00,.9448967E+00/
      data (xp(i),i= 921, 960)/
     $.9465965E+00,.9482703E+00,.9499181E+00,.9515400E+00,.9531358E+00,
     $.9547056E+00,.9562492E+00,.9577666E+00,.9592578E+00,.9607228E+00,
     $.9621615E+00,.9635738E+00,.9649597E+00,.9663193E+00,.9676524E+00,
     $.9689590E+00,.9702391E+00,.9714927E+00,.9727196E+00,.9739199E+00,
     $.9750936E+00,.9762406E+00,.9773609E+00,.9784544E+00,.9795211E+00,
     $.9805610E+00,.9815741E+00,.9825603E+00,.9835197E+00,.9844521E+00,
     $.9853575E+00,.9862360E+00,.9870876E+00,.9879120E+00,.9887095E+00,
     $.9894799E+00,.9902232E+00,.9909394E+00,.9916285E+00,.9922904E+00/
      data (xp(i),i= 961,1000)/
     $.9929252E+00,.9935328E+00,.9941132E+00,.9946664E+00,.9951924E+00,
     $.9956911E+00,.9961626E+00,.9966068E+00,.9970238E+00,.9974135E+00,
     $.9977758E+00,.9981109E+00,.9984186E+00,.9986990E+00,.9989521E+00,
     $.9991778E+00,.9993762E+00,.9995472E+00,.9996909E+00,.9998072E+00,
     $.9998962E+00,.9999577E+00,.9999920E+00,.1962267E-02,.5886772E-02,
     $.9811186E-02,.1373545E-01,.1765950E-01,.2158328E-01,.2550673E-01,
     $.2942978E-01,.3335238E-01,.3727447E-01,.4119598E-01,.4511686E-01,
     $.4903704E-01,.5295647E-01,.5687508E-01,.6079282E-01,.6470962E-01/
      data (xp(i),i=1001,1040)/
     $.6862542E-01,.7254017E-01,.7645380E-01,.8036625E-01,.8427746E-01,
     $.8818738E-01,.9209594E-01,.9600308E-01,.9990874E-01,.1038129E+00,
     $.1077154E+00,.1116162E+00,.1155154E+00,.1194128E+00,.1233083E+00,
     $.1272019E+00,.1310936E+00,.1349832E+00,.1388708E+00,.1427562E+00,
     $.1466395E+00,.1505204E+00,.1543991E+00,.1582754E+00,.1621492E+00,
     $.1660205E+00,.1698893E+00,.1737555E+00,.1776190E+00,.1814798E+00,
     $.1853377E+00,.1891928E+00,.1930450E+00,.1968942E+00,.2007404E+00,
     $.2045835E+00,.2084235E+00,.2122602E+00,.2160937E+00,.2199238E+00/
      data (xp(i),i=1041,1080)/
     $.2237505E+00,.2275738E+00,.2313936E+00,.2352099E+00,.2390225E+00,
     $.2428314E+00,.2466366E+00,.2504380E+00,.2542355E+00,.2580292E+00,
     $.2618188E+00,.2656044E+00,.2693859E+00,.2731633E+00,.2769365E+00,
     $.2807054E+00,.2844699E+00,.2882301E+00,.2919859E+00,.2957372E+00,
     $.2994839E+00,.3032259E+00,.3069634E+00,.3106961E+00,.3144240E+00,
     $.3181470E+00,.3218652E+00,.3255784E+00,.3292866E+00,.3329897E+00,
     $.3366877E+00,.3403805E+00,.3440681E+00,.3477503E+00,.3514272E+00,
     $.3550987E+00,.3587648E+00,.3624253E+00,.3660802E+00,.3697295E+00/
      data (xp(i),i=1081,1120)/
     $.3733731E+00,.3770109E+00,.3806429E+00,.3842691E+00,.3878893E+00,
     $.3915036E+00,.3951118E+00,.3987140E+00,.4023100E+00,.4058998E+00,
     $.4094834E+00,.4130607E+00,.4166316E+00,.4201960E+00,.4237541E+00,
     $.4273055E+00,.4308504E+00,.4343887E+00,.4379203E+00,.4414451E+00,
     $.4449632E+00,.4484743E+00,.4519786E+00,.4554759E+00,.4589662E+00,
     $.4624494E+00,.4659255E+00,.4693945E+00,.4728562E+00,.4763106E+00,
     $.4797577E+00,.4831973E+00,.4866296E+00,.4900543E+00,.4934715E+00,
     $.4968812E+00,.5002831E+00,.5036774E+00,.5070638E+00,.5104425E+00/
      data (xp(i),i=1121,1160)/
     $.5138133E+00,.5171762E+00,.5205312E+00,.5238781E+00,.5272169E+00,
     $.5305477E+00,.5338702E+00,.5371846E+00,.5404906E+00,.5437884E+00,
     $.5470777E+00,.5503587E+00,.5536311E+00,.5568951E+00,.5601504E+00,
     $.5633972E+00,.5666352E+00,.5698645E+00,.5730851E+00,.5762968E+00,
     $.5794996E+00,.5826936E+00,.5858785E+00,.5890544E+00,.5922213E+00,
     $.5953790E+00,.5985276E+00,.6016669E+00,.6047970E+00,.6079177E+00,
     $.6110291E+00,.6141311E+00,.6172236E+00,.6203066E+00,.6233801E+00,
     $.6264440E+00,.6294982E+00,.6325427E+00,.6355775E+00,.6386024E+00/
      data (xp(i),i=1161,1200)/
     $.6416176E+00,.6446229E+00,.6476182E+00,.6506036E+00,.6535789E+00,
     $.6565442E+00,.6594994E+00,.6624444E+00,.6653792E+00,.6683037E+00,
     $.6712180E+00,.6741219E+00,.6770155E+00,.6798986E+00,.6827712E+00,
     $.6856333E+00,.6884849E+00,.6913259E+00,.6941562E+00,.6969758E+00,
     $.6997847E+00,.7025828E+00,.7053701E+00,.7081465E+00,.7109120E+00,
     $.7136666E+00,.7164102E+00,.7191427E+00,.7218642E+00,.7245746E+00,
     $.7272737E+00,.7299617E+00,.7326385E+00,.7353039E+00,.7379581E+00,
     $.7406008E+00,.7432322E+00,.7458521E+00,.7484606E+00,.7510575E+00/
      data (xp(i),i=1201,1240)/
     $.7536428E+00,.7562165E+00,.7587786E+00,.7613290E+00,.7638676E+00,
     $.7663945E+00,.7689096E+00,.7714129E+00,.7739043E+00,.7763837E+00,
     $.7788512E+00,.7813067E+00,.7837502E+00,.7861816E+00,.7886009E+00,
     $.7910080E+00,.7934030E+00,.7957857E+00,.7981562E+00,.8005144E+00,
     $.8028602E+00,.8051937E+00,.8075148E+00,.8098234E+00,.8121196E+00,
     $.8144033E+00,.8166744E+00,.8189330E+00,.8211789E+00,.8234122E+00,
     $.8256328E+00,.8278407E+00,.8300358E+00,.8322182E+00,.8343877E+00,
     $.8365444E+00,.8386882E+00,.8408191E+00,.8429370E+00,.8450420E+00/
      data (xp(i),i=1241,1280)/
     $.8471339E+00,.8492128E+00,.8512786E+00,.8533313E+00,.8553709E+00,
     $.8573972E+00,.8594104E+00,.8614104E+00,.8633970E+00,.8653704E+00,
     $.8673304E+00,.8692771E+00,.8712104E+00,.8731303E+00,.8750367E+00,
     $.8769297E+00,.8788091E+00,.8806750E+00,.8825274E+00,.8843662E+00,
     $.8861913E+00,.8880028E+00,.8898006E+00,.8915847E+00,.8933550E+00,
     $.8951117E+00,.8968545E+00,.8985835E+00,.9002987E+00,.9020000E+00,
     $.9036874E+00,.9053609E+00,.9070204E+00,.9086660E+00,.9102976E+00,
     $.9119152E+00,.9135187E+00,.9151081E+00,.9166835E+00,.9182447E+00/
      data (xp(i),i=1281,1320)/
     $.9197918E+00,.9213247E+00,.9228435E+00,.9243480E+00,.9258383E+00,
     $.9273143E+00,.9287760E+00,.9302235E+00,.9316566E+00,.9330754E+00,
     $.9344797E+00,.9358697E+00,.9372453E+00,.9386065E+00,.9399532E+00,
     $.9412854E+00,.9426031E+00,.9439063E+00,.9451950E+00,.9464691E+00,
     $.9477286E+00,.9489735E+00,.9502038E+00,.9514195E+00,.9526205E+00,
     $.9538069E+00,.9549785E+00,.9561355E+00,.9572777E+00,.9584052E+00,
     $.9595179E+00,.9606159E+00,.9616990E+00,.9627673E+00,.9638208E+00,
     $.9648595E+00,.9658833E+00,.9668922E+00,.9678863E+00,.9688654E+00/
      data (xp(i),i=1321,1360)/
     $.9698296E+00,.9707788E+00,.9717132E+00,.9726325E+00,.9735369E+00,
     $.9744262E+00,.9753006E+00,.9761600E+00,.9770043E+00,.9778335E+00,
     $.9786477E+00,.9794468E+00,.9802309E+00,.9809998E+00,.9817537E+00,
     $.9824924E+00,.9832160E+00,.9839244E+00,.9846177E+00,.9852958E+00,
     $.9859587E+00,.9866065E+00,.9872390E+00,.9878564E+00,.9884585E+00,
     $.9890455E+00,.9896171E+00,.9901736E+00,.9907148E+00,.9912407E+00,
     $.9917514E+00,.9922468E+00,.9927269E+00,.9931917E+00,.9936412E+00,
     $.9940754E+00,.9944943E+00,.9948979E+00,.9952862E+00,.9956591E+00/
      data (xp(i),i=1361,1400)/
     $.9960167E+00,.9963590E+00,.9966859E+00,.9969974E+00,.9972936E+00,
     $.9975745E+00,.9978400E+00,.9980901E+00,.9983248E+00,.9985442E+00,
     $.9987482E+00,.9989368E+00,.9991100E+00,.9992678E+00,.9994103E+00,
     $.9995373E+00,.9996489E+00,.9997452E+00,.9998261E+00,.9998915E+00,
     $.9999416E+00,.9999762E+00,.9999955E+00,.1570010E-02,.4710016E-02,
     $.7849975E-02,.1098986E-01,.1412963E-01,.1726926E-01,.2040873E-01,
     $.2354799E-01,.2668702E-01,.2982579E-01,.3296426E-01,.3610241E-01,
     $.3924020E-01,.4237761E-01,.4551459E-01,.4865113E-01,.5178719E-01/
      data (xp(i),i=1401,1440)/
     $.5492274E-01,.5805775E-01,.6119218E-01,.6432601E-01,.6745921E-01,
     $.7059174E-01,.7372358E-01,.7685468E-01,.7998504E-01,.8311460E-01,
     $.8624334E-01,.8937123E-01,.9249824E-01,.9562434E-01,.9874950E-01,
     $.1018737E+00,.1049969E+00,.1081190E+00,.1112401E+00,.1143601E+00,
     $.1174789E+00,.1205966E+00,.1237131E+00,.1268284E+00,.1299424E+00,
     $.1330552E+00,.1361666E+00,.1392767E+00,.1423855E+00,.1454928E+00,
     $.1485987E+00,.1517031E+00,.1548060E+00,.1579074E+00,.1610073E+00,
     $.1641055E+00,.1672022E+00,.1702971E+00,.1733905E+00,.1764821E+00/
      data (xp(i),i=1441,1480)/
     $.1795719E+00,.1826600E+00,.1857463E+00,.1888308E+00,.1919134E+00,
     $.1949941E+00,.1980728E+00,.2011497E+00,.2042245E+00,.2072973E+00,
     $.2103681E+00,.2134368E+00,.2165035E+00,.2195679E+00,.2226302E+00,
     $.2256904E+00,.2287482E+00,.2318039E+00,.2348572E+00,.2379083E+00,
     $.2409569E+00,.2440033E+00,.2470472E+00,.2500886E+00,.2531276E+00,
     $.2561641E+00,.2591981E+00,.2622295E+00,.2652584E+00,.2682846E+00,
     $.2713082E+00,.2743291E+00,.2773473E+00,.2803628E+00,.2833755E+00,
     $.2863854E+00,.2893925E+00,.2923967E+00,.2953980E+00,.2983965E+00/
      data (xp(i),i=1481,1520)/
     $.3013920E+00,.3043845E+00,.3073740E+00,.3103605E+00,.3133439E+00,
     $.3163243E+00,.3193015E+00,.3222756E+00,.3252465E+00,.3282141E+00,
     $.3311786E+00,.3341398E+00,.3370977E+00,.3400522E+00,.3430035E+00,
     $.3459513E+00,.3488957E+00,.3518367E+00,.3547742E+00,.3577082E+00,
     $.3606387E+00,.3635657E+00,.3664890E+00,.3694087E+00,.3723248E+00,
     $.3752373E+00,.3781460E+00,.3810510E+00,.3839522E+00,.3868497E+00,
     $.3897433E+00,.3926331E+00,.3955190E+00,.3984010E+00,.4012791E+00,
     $.4041533E+00,.4070234E+00,.4098896E+00,.4127517E+00,.4156097E+00/
      data (xp(i),i=1521,1560)/
     $.4184636E+00,.4213134E+00,.4241591E+00,.4270006E+00,.4298378E+00,
     $.4326708E+00,.4354996E+00,.4383241E+00,.4411442E+00,.4439600E+00,
     $.4467714E+00,.4495784E+00,.4523810E+00,.4551791E+00,.4579727E+00,
     $.4607618E+00,.4635464E+00,.4663264E+00,.4691018E+00,.4718726E+00,
     $.4746387E+00,.4774001E+00,.4801569E+00,.4829089E+00,.4856561E+00,
     $.4883986E+00,.4911362E+00,.4938690E+00,.4965969E+00,.4993199E+00,
     $.5020381E+00,.5047512E+00,.5074594E+00,.5101626E+00,.5128607E+00,
     $.5155538E+00,.5182418E+00,.5209247E+00,.5236025E+00,.5262750E+00/
      data (xp(i),i=1561,1600)/
     $.5289425E+00,.5316046E+00,.5342616E+00,.5369133E+00,.5395597E+00,
     $.5422007E+00,.5448365E+00,.5474668E+00,.5500918E+00,.5527113E+00,
     $.5553254E+00,.5579340E+00,.5605371E+00,.5631347E+00,.5657267E+00,
     $.5683131E+00,.5708940E+00,.5734692E+00,.5760387E+00,.5786026E+00,
     $.5811608E+00,.5837132E+00,.5862599E+00,.5888008E+00,.5913360E+00,
     $.5938652E+00,.5963886E+00,.5989062E+00,.6014178E+00,.6039235E+00,
     $.6064233E+00,.6089170E+00,.6114048E+00,.6138865E+00,.6163622E+00,
     $.6188318E+00,.6212953E+00,.6237527E+00,.6262040E+00,.6286490E+00/
      data (xp(i),i=1601,1640)/
     $.6310879E+00,.6335205E+00,.6359469E+00,.6383670E+00,.6407808E+00,
     $.6431883E+00,.6455895E+00,.6479843E+00,.6503727E+00,.6527547E+00,
     $.6551303E+00,.6574994E+00,.6598620E+00,.6622181E+00,.6645677E+00,
     $.6669107E+00,.6692472E+00,.6715770E+00,.6739003E+00,.6762169E+00,
     $.6785268E+00,.6808300E+00,.6831266E+00,.6854163E+00,.6876994E+00,
     $.6899756E+00,.6922451E+00,.6945077E+00,.6967635E+00,.6990124E+00,
     $.7012544E+00,.7034895E+00,.7057176E+00,.7079388E+00,.7101530E+00,
     $.7123603E+00,.7145605E+00,.7167536E+00,.7189397E+00,.7211187E+00/
      data (xp(i),i=1641,1680)/
     $.7232906E+00,.7254553E+00,.7276129E+00,.7297634E+00,.7319066E+00,
     $.7340426E+00,.7361714E+00,.7382929E+00,.7404071E+00,.7425141E+00,
     $.7446137E+00,.7467060E+00,.7487909E+00,.7508684E+00,.7529385E+00,
     $.7550013E+00,.7570565E+00,.7591043E+00,.7611446E+00,.7631774E+00,
     $.7652027E+00,.7672204E+00,.7692306E+00,.7712332E+00,.7732282E+00,
     $.7752155E+00,.7771953E+00,.7791673E+00,.7811317E+00,.7830884E+00,
     $.7850373E+00,.7869785E+00,.7889120E+00,.7908376E+00,.7927555E+00,
     $.7946656E+00,.7965678E+00,.7984622E+00,.8003486E+00,.8022273E+00/
      data (xp(i),i=1681,1720)/
     $.8040979E+00,.8059607E+00,.8078155E+00,.8096624E+00,.8115013E+00,
     $.8133321E+00,.8151550E+00,.8169698E+00,.8187765E+00,.8205752E+00,
     $.8223658E+00,.8241483E+00,.8259227E+00,.8276889E+00,.8294470E+00,
     $.8311968E+00,.8329385E+00,.8346720E+00,.8363972E+00,.8381142E+00,
     $.8398229E+00,.8415234E+00,.8432156E+00,.8448994E+00,.8465749E+00,
     $.8482421E+00,.8499009E+00,.8515513E+00,.8531933E+00,.8548269E+00,
     $.8564521E+00,.8580688E+00,.8596771E+00,.8612769E+00,.8628682E+00,
     $.8644510E+00,.8660253E+00,.8675910E+00,.8691482E+00,.8706968E+00/
      data (xp(i),i=1721,1760)/
     $.8722369E+00,.8737683E+00,.8752911E+00,.8768053E+00,.8783108E+00,
     $.8798077E+00,.8812959E+00,.8827754E+00,.8842463E+00,.8857083E+00,
     $.8871617E+00,.8886063E+00,.8900422E+00,.8914692E+00,.8928875E+00,
     $.8942970E+00,.8956977E+00,.8970895E+00,.8984725E+00,.8998466E+00,
     $.9012119E+00,.9025683E+00,.9039157E+00,.9052543E+00,.9065839E+00,
     $.9079046E+00,.9092164E+00,.9105192E+00,.9118130E+00,.9130978E+00,
     $.9143736E+00,.9156404E+00,.9168982E+00,.9181469E+00,.9193866E+00,
     $.9206172E+00,.9218387E+00,.9230511E+00,.9242545E+00,.9254487E+00/
      data (xp(i),i=1761,1800)/
     $.9266338E+00,.9278098E+00,.9289766E+00,.9301343E+00,.9312828E+00,
     $.9324221E+00,.9335522E+00,.9346731E+00,.9357848E+00,.9368872E+00,
     $.9379805E+00,.9390645E+00,.9401392E+00,.9412046E+00,.9422608E+00,
     $.9433077E+00,.9443453E+00,.9453735E+00,.9463925E+00,.9474021E+00,
     $.9484024E+00,.9493933E+00,.9503749E+00,.9513471E+00,.9523099E+00,
     $.9532633E+00,.9542073E+00,.9551420E+00,.9560672E+00,.9569829E+00,
     $.9578893E+00,.9587862E+00,.9596736E+00,.9605516E+00,.9614201E+00,
     $.9622791E+00,.9631287E+00,.9639687E+00,.9647992E+00,.9656203E+00/
      data (xp(i),i=1801,1840)/
     $.9664318E+00,.9672338E+00,.9680262E+00,.9688091E+00,.9695824E+00,
     $.9703462E+00,.9711004E+00,.9718451E+00,.9725801E+00,.9733056E+00,
     $.9740215E+00,.9747278E+00,.9754244E+00,.9761115E+00,.9767889E+00,
     $.9774567E+00,.9781148E+00,.9787633E+00,.9794022E+00,.9800314E+00,
     $.9806509E+00,.9812608E+00,.9818610E+00,.9824515E+00,.9830323E+00,
     $.9836035E+00,.9841649E+00,.9847166E+00,.9852586E+00,.9857909E+00,
     $.9863135E+00,.9868264E+00,.9873295E+00,.9878229E+00,.9883066E+00,
     $.9887805E+00,.9892447E+00,.9896991E+00,.9901437E+00,.9905786E+00/
      data (xp(i),i=1841,1880)/
     $.9910037E+00,.9914191E+00,.9918247E+00,.9922205E+00,.9926065E+00,
     $.9929827E+00,.9933492E+00,.9937058E+00,.9940527E+00,.9943897E+00,
     $.9947169E+00,.9950344E+00,.9953420E+00,.9956398E+00,.9959278E+00,
     $.9962060E+00,.9964743E+00,.9967328E+00,.9969815E+00,.9972204E+00,
     $.9974494E+00,.9976686E+00,.9978780E+00,.9980775E+00,.9982672E+00,
     $.9984471E+00,.9986171E+00,.9987772E+00,.9989275E+00,.9990680E+00,
     $.9991986E+00,.9993193E+00,.9994302E+00,.9995313E+00,.9996225E+00,
     $.9997038E+00,.9997753E+00,.9998369E+00,.9998886E+00,.9999306E+00/
      data (xp(i),i=1881,1883)/
     $.9999626E+00,.9999848E+00,.9999971E+00/
      data (wp(i),i=   1,  40)/
     $.4679139E+00,.3607616E+00,.1713245E+00,.2955242E+00,.2692667E+00,
     $.2190864E+00,.1494513E+00,.6667134E-01,.1527534E+00,.1491730E+00,
     $.1420961E+00,.1316886E+00,.1181945E+00,.1019301E+00,.8327674E-01,
     $.6267205E-01,.4060143E-01,.1761401E-01,.7750595E-01,.7703982E-01,
     $.7611036E-01,.7472317E-01,.7288658E-01,.7061165E-01,.6791205E-01,
     $.6480401E-01,.6130624E-01,.5743977E-01,.5322785E-01,.4869581E-01,
     $.4387091E-01,.3878217E-01,.3346020E-01,.2793701E-01,.2224585E-01,
     $.1642106E-01,.1049828E-01,.4521277E-02,.5190788E-01,.5176794E-01/
      data (wp(i),i=  41,  80)/
     $.5148845E-01,.5107016E-01,.5051418E-01,.4982204E-01,.4899558E-01,
     $.4803703E-01,.4694899E-01,.4573438E-01,.4439648E-01,.4293889E-01,
     $.4136555E-01,.3968070E-01,.3788887E-01,.3599490E-01,.3400389E-01,
     $.3192122E-01,.2975249E-01,.2750356E-01,.2518048E-01,.2278952E-01,
     $.2033712E-01,.1782990E-01,.1527462E-01,.1267817E-01,.1004756E-01,
     $.7389931E-02,.4712730E-02,.2026812E-02,.3901781E-01,.3895840E-01,
     $.3883965E-01,.3866176E-01,.3842499E-01,.3812971E-01,.3777636E-01,
     $.3736549E-01,.3689771E-01,.3637375E-01,.3579439E-01,.3516053E-01/
      data (wp(i),i=  81, 120)/
     $.3447312E-01,.3373321E-01,.3294194E-01,.3210050E-01,.3121017E-01,
     $.3027232E-01,.2928837E-01,.2825982E-01,.2718823E-01,.2607524E-01,
     $.2492254E-01,.2373188E-01,.2250509E-01,.2124403E-01,.1995061E-01,
     $.1862681E-01,.1727465E-01,.1589618E-01,.1449351E-01,.1306876E-01,
     $.1162411E-01,.1016177E-01,.8683945E-02,.7192905E-02,.5690922E-02,
     $.4180313E-02,.2663534E-02,.1144950E-02,.3125542E-01,.3122488E-01,
     $.3116384E-01,.3107234E-01,.3095048E-01,.3079838E-01,.3061619E-01,
     $.3040408E-01,.3016227E-01,.2989098E-01,.2959049E-01,.2926108E-01/
      data (wp(i),i= 121, 160)/
     $.2890309E-01,.2851685E-01,.2810276E-01,.2766120E-01,.2719261E-01,
     $.2669746E-01,.2617622E-01,.2562940E-01,.2505754E-01,.2446120E-01,
     $.2384096E-01,.2319742E-01,.2253122E-01,.2184300E-01,.2113344E-01,
     $.2040323E-01,.1965309E-01,.1888374E-01,.1809594E-01,.1729046E-01,
     $.1646809E-01,.1562962E-01,.1477588E-01,.1390771E-01,.1302595E-01,
     $.1213146E-01,.1122511E-01,.1030780E-01,.9380420E-02,.8443871E-02,
     $.7499073E-02,.6546948E-02,.5588428E-02,.4624450E-02,.3655961E-02,
     $.2683925E-02,.1709393E-02,.7346345E-03,.2087312E-01,.2086402E-01/
      data (wp(i),i= 161, 200)/
     $.2084584E-01,.2081857E-01,.2078223E-01,.2073683E-01,.2068240E-01,
     $.2061896E-01,.2054653E-01,.2046515E-01,.2037486E-01,.2027568E-01,
     $.2016767E-01,.2005088E-01,.1992534E-01,.1979113E-01,.1964829E-01,
     $.1949689E-01,.1933700E-01,.1916867E-01,.1899200E-01,.1880705E-01,
     $.1861391E-01,.1841266E-01,.1820338E-01,.1798617E-01,.1776113E-01,
     $.1752835E-01,.1728792E-01,.1703997E-01,.1678459E-01,.1652190E-01,
     $.1625201E-01,.1597504E-01,.1569110E-01,.1540033E-01,.1510285E-01,
     $.1479879E-01,.1448828E-01,.1417146E-01,.1384846E-01,.1351943E-01/
      data (wp(i),i= 201, 240)/
     $.1318451E-01,.1284384E-01,.1249758E-01,.1214587E-01,.1178887E-01,
     $.1142673E-01,.1105962E-01,.1068768E-01,.1031109E-01,.9930004E-02,
     $.9544593E-02,.9155022E-02,.8761463E-02,.8364086E-02,.7963064E-02,
     $.7558573E-02,.7150788E-02,.6739888E-02,.6326051E-02,.5909457E-02,
     $.5490289E-02,.5068728E-02,.4644959E-02,.4219166E-02,.3791535E-02,
     $.3362252E-02,.2931504E-02,.2499479E-02,.2066366E-02,.1632357E-02,
     $.1197647E-02,.7624721E-03,.3276087E-03,.1566826E-01,.1566442E-01,
     $.1565672E-01,.1564519E-01,.1562981E-01,.1561059E-01,.1558755E-01/
      data (wp(i),i= 241, 280)/
     $.1556067E-01,.1552998E-01,.1549547E-01,.1545716E-01,.1541506E-01,
     $.1536917E-01,.1531950E-01,.1526608E-01,.1520891E-01,.1514800E-01,
     $.1508338E-01,.1501505E-01,.1494303E-01,.1486735E-01,.1478802E-01,
     $.1470505E-01,.1461848E-01,.1452832E-01,.1443459E-01,.1433731E-01,
     $.1423652E-01,.1413223E-01,.1402447E-01,.1391327E-01,.1379866E-01,
     $.1368065E-01,.1355929E-01,.1343460E-01,.1330661E-01,.1317535E-01,
     $.1304086E-01,.1290316E-01,.1276230E-01,.1261831E-01,.1247122E-01,
     $.1232106E-01,.1216788E-01,.1201172E-01,.1185260E-01,.1169058E-01/
      data (wp(i),i= 281, 320)/
     $.1152568E-01,.1135796E-01,.1118744E-01,.1101418E-01,.1083822E-01,
     $.1065959E-01,.1047835E-01,.1029454E-01,.1010820E-01,.9919373E-02,
     $.9728115E-02,.9534468E-02,.9338480E-02,.9140200E-02,.8939676E-02,
     $.8736957E-02,.8532093E-02,.8325134E-02,.8116132E-02,.7905137E-02,
     $.7692201E-02,.7477377E-02,.7260717E-02,.7042274E-02,.6822103E-02,
     $.6600256E-02,.6376790E-02,.6151757E-02,.5925215E-02,.5697218E-02,
     $.5467822E-02,.5237083E-02,.5005059E-02,.4771806E-02,.4537382E-02,
     $.4301844E-02,.4065249E-02,.3827657E-02,.3589125E-02,.3349711E-02/
      data (wp(i),i= 321, 360)/
     $.3109476E-02,.2868477E-02,.2626773E-02,.2384425E-02,.2141492E-02,
     $.1898033E-02,.1654108E-02,.1409777E-02,.1165101E-02,.9201405E-03,
     $.6749606E-03,.4296466E-03,.1845901E-03,.1045439E-01,.1045325E-01,
     $.1045097E-01,.1044754E-01,.1044297E-01,.1043726E-01,.1043041E-01,
     $.1042242E-01,.1041329E-01,.1040302E-01,.1039161E-01,.1037907E-01,
     $.1036539E-01,.1035058E-01,.1033464E-01,.1031758E-01,.1029938E-01,
     $.1028006E-01,.1025961E-01,.1023804E-01,.1021535E-01,.1019155E-01,
     $.1016663E-01,.1014060E-01,.1011347E-01,.1008523E-01,.1005588E-01/
      data (wp(i),i= 361, 400)/
     $.1002544E-01,.9993899E-02,.9961267E-02,.9927547E-02,.9892741E-02,
     $.9856855E-02,.9819891E-02,.9781854E-02,.9742747E-02,.9702576E-02,
     $.9661345E-02,.9619057E-02,.9575718E-02,.9531333E-02,.9485905E-02,
     $.9439441E-02,.9391946E-02,.9343424E-02,.9293880E-02,.9243321E-02,
     $.9191751E-02,.9139177E-02,.9085604E-02,.9031038E-02,.8975485E-02,
     $.8918951E-02,.8861442E-02,.8802965E-02,.8743525E-02,.8683130E-02,
     $.8621786E-02,.8559499E-02,.8496277E-02,.8432127E-02,.8367054E-02,
     $.8301068E-02,.8234174E-02,.8166380E-02,.8097693E-02,.8028121E-02/
      data (wp(i),i= 401, 440)/
     $.7957672E-02,.7886353E-02,.7814173E-02,.7741138E-02,.7667257E-02,
     $.7592538E-02,.7516989E-02,.7440619E-02,.7363435E-02,.7285447E-02,
     $.7206662E-02,.7127090E-02,.7046739E-02,.6965617E-02,.6883734E-02,
     $.6801099E-02,.6717721E-02,.6633608E-02,.6548770E-02,.6463217E-02,
     $.6376957E-02,.6290000E-02,.6202356E-02,.6114033E-02,.6025043E-02,
     $.5935394E-02,.5845096E-02,.5754159E-02,.5662594E-02,.5570409E-02,
     $.5477616E-02,.5384224E-02,.5290244E-02,.5195685E-02,.5100559E-02,
     $.5004875E-02,.4908644E-02,.4811876E-02,.4714583E-02,.4616774E-02/
      data (wp(i),i= 441, 480)/
     $.4518461E-02,.4419654E-02,.4320364E-02,.4220601E-02,.4120378E-02,
     $.4019704E-02,.3918590E-02,.3817049E-02,.3715090E-02,.3612725E-02,
     $.3509965E-02,.3406822E-02,.3303306E-02,.3199429E-02,.3095203E-02,
     $.2990638E-02,.2885746E-02,.2780539E-02,.2675029E-02,.2569225E-02,
     $.2463141E-02,.2356788E-02,.2250177E-02,.2143320E-02,.2036229E-02,
     $.1928915E-02,.1821391E-02,.1713667E-02,.1605756E-02,.1497670E-02,
     $.1389420E-02,.1281018E-02,.1172476E-02,.1063806E-02,.9550200E-03,
     $.8461294E-03,.7371464E-03,.6280830E-03,.5189512E-03,.4097636E-03/
      data (wp(i),i= 481, 520)/
     $.3005340E-03,.1912855E-03,.8217779E-04,.7844110E-02,.7843627E-02,
     $.7842662E-02,.7841214E-02,.7839284E-02,.7836871E-02,.7833976E-02,
     $.7830599E-02,.7826741E-02,.7822400E-02,.7817579E-02,.7812276E-02,
     $.7806493E-02,.7800229E-02,.7793485E-02,.7786262E-02,.7778560E-02,
     $.7770379E-02,.7761720E-02,.7752583E-02,.7742970E-02,.7732880E-02,
     $.7722314E-02,.7711273E-02,.7699757E-02,.7687768E-02,.7675306E-02,
     $.7662371E-02,.7648965E-02,.7635088E-02,.7620742E-02,.7605926E-02,
     $.7590643E-02,.7574892E-02,.7558676E-02,.7541994E-02,.7524848E-02/
      data (wp(i),i= 521, 560)/
     $.7507240E-02,.7489169E-02,.7470638E-02,.7451646E-02,.7432197E-02,
     $.7412290E-02,.7391927E-02,.7371109E-02,.7349838E-02,.7328114E-02,
     $.7305939E-02,.7283315E-02,.7260243E-02,.7236724E-02,.7212760E-02,
     $.7188352E-02,.7163501E-02,.7138210E-02,.7112480E-02,.7086312E-02,
     $.7059708E-02,.7032669E-02,.7005198E-02,.6977296E-02,.6948964E-02,
     $.6920205E-02,.6891020E-02,.6861412E-02,.6831380E-02,.6800929E-02,
     $.6770059E-02,.6738773E-02,.6707072E-02,.6674958E-02,.6642433E-02,
     $.6609500E-02,.6576160E-02,.6542416E-02,.6508269E-02,.6473721E-02/
      data (wp(i),i= 561, 600)/
     $.6438775E-02,.6403433E-02,.6367697E-02,.6331569E-02,.6295052E-02,
     $.6258147E-02,.6220857E-02,.6183184E-02,.6145131E-02,.6106700E-02,
     $.6067893E-02,.6028713E-02,.5989161E-02,.5949241E-02,.5908956E-02,
     $.5868306E-02,.5827296E-02,.5785926E-02,.5744201E-02,.5702123E-02,
     $.5659693E-02,.5616916E-02,.5573792E-02,.5530326E-02,.5486520E-02,
     $.5442376E-02,.5397897E-02,.5353085E-02,.5307945E-02,.5262478E-02,
     $.5216687E-02,.5170575E-02,.5124145E-02,.5077400E-02,.5030342E-02,
     $.4982975E-02,.4935301E-02,.4887323E-02,.4839045E-02,.4790469E-02/
      data (wp(i),i= 601, 640)/
     $.4741598E-02,.4692436E-02,.4642984E-02,.4593247E-02,.4543228E-02,
     $.4492929E-02,.4442353E-02,.4391504E-02,.4340385E-02,.4288999E-02,
     $.4237349E-02,.4185438E-02,.4133270E-02,.4080847E-02,.4028173E-02,
     $.3975251E-02,.3922085E-02,.3868678E-02,.3815032E-02,.3761152E-02,
     $.3707040E-02,.3652700E-02,.3598135E-02,.3543349E-02,.3488345E-02,
     $.3433126E-02,.3377697E-02,.3322059E-02,.3266217E-02,.3210173E-02,
     $.3153933E-02,.3097498E-02,.3040873E-02,.2984060E-02,.2927064E-02,
     $.2869888E-02,.2812535E-02,.2755009E-02,.2697314E-02,.2639453E-02/
      data (wp(i),i= 641, 680)/
     $.2581429E-02,.2523247E-02,.2464909E-02,.2406419E-02,.2347782E-02,
     $.2289000E-02,.2230077E-02,.2171017E-02,.2111823E-02,.2052500E-02,
     $.1993050E-02,.1933477E-02,.1873786E-02,.1813979E-02,.1754061E-02,
     $.1694034E-02,.1633904E-02,.1573673E-02,.1513345E-02,.1452924E-02,
     $.1392413E-02,.1331817E-02,.1271139E-02,.1210383E-02,.1149552E-02,
     $.1088651E-02,.1027682E-02,.9666507E-03,.9055595E-03,.8444126E-03,
     $.7832138E-03,.7219667E-03,.6606753E-03,.5993432E-03,.5379742E-03,
     $.4765722E-03,.4151409E-03,.3536841E-03,.2922057E-03,.2307099E-03/
      data (wp(i),i= 681, 720)/
     $.1692014E-03,.1076904E-03,.4626372E-04,.5231608E-02,.5231465E-02,
     $.5231179E-02,.5230749E-02,.5230177E-02,.5229461E-02,.5228602E-02,
     $.5227600E-02,.5226454E-02,.5225166E-02,.5223735E-02,.5222161E-02,
     $.5220444E-02,.5218584E-02,.5216581E-02,.5214435E-02,.5212147E-02,
     $.5209716E-02,.5207142E-02,.5204426E-02,.5201567E-02,.5198567E-02,
     $.5195423E-02,.5192138E-02,.5188710E-02,.5185141E-02,.5181429E-02,
     $.5177576E-02,.5173581E-02,.5169445E-02,.5165167E-02,.5160747E-02,
     $.5156186E-02,.5151485E-02,.5146642E-02,.5141658E-02,.5136534E-02/
      data (wp(i),i= 721, 760)/
     $.5131269E-02,.5125863E-02,.5120318E-02,.5114632E-02,.5108806E-02,
     $.5102840E-02,.5096735E-02,.5090490E-02,.5084106E-02,.5077583E-02,
     $.5070920E-02,.5064119E-02,.5057180E-02,.5050102E-02,.5042885E-02,
     $.5035531E-02,.5028039E-02,.5020409E-02,.5012642E-02,.5004738E-02,
     $.4996696E-02,.4988518E-02,.4980203E-02,.4971752E-02,.4963165E-02,
     $.4954443E-02,.4945584E-02,.4936590E-02,.4927461E-02,.4918197E-02,
     $.4908799E-02,.4899266E-02,.4889599E-02,.4879799E-02,.4869864E-02,
     $.4859797E-02,.4849596E-02,.4839263E-02,.4828797E-02,.4818199E-02/
      data (wp(i),i= 761, 800)/
     $.4807470E-02,.4796608E-02,.4785616E-02,.4774492E-02,.4763238E-02,
     $.4751853E-02,.4740338E-02,.4728694E-02,.4716920E-02,.4705017E-02,
     $.4692985E-02,.4680825E-02,.4668537E-02,.4656121E-02,.4643577E-02,
     $.4630907E-02,.4618109E-02,.4605185E-02,.4592136E-02,.4578960E-02,
     $.4565659E-02,.4552233E-02,.4538683E-02,.4525008E-02,.4511210E-02,
     $.4497288E-02,.4483243E-02,.4469075E-02,.4454785E-02,.4440373E-02,
     $.4425840E-02,.4411185E-02,.4396410E-02,.4381514E-02,.4366498E-02,
     $.4351363E-02,.4336109E-02,.4320736E-02,.4305245E-02,.4289636E-02/
      data (wp(i),i= 801, 840)/
     $.4273910E-02,.4258066E-02,.4242106E-02,.4226030E-02,.4209839E-02,
     $.4193532E-02,.4177110E-02,.4160574E-02,.4143924E-02,.4127161E-02,
     $.4110284E-02,.4093296E-02,.4076195E-02,.4058982E-02,.4041659E-02,
     $.4024225E-02,.4006681E-02,.3989027E-02,.3971264E-02,.3953392E-02,
     $.3935412E-02,.3917324E-02,.3899129E-02,.3880828E-02,.3862420E-02,
     $.3843906E-02,.3825288E-02,.3806564E-02,.3787737E-02,.3768805E-02,
     $.3749771E-02,.3730634E-02,.3711395E-02,.3692054E-02,.3672612E-02,
     $.3653070E-02,.3633427E-02,.3613685E-02,.3593845E-02,.3573906E-02/
      data (wp(i),i= 841, 880)/
     $.3553869E-02,.3533735E-02,.3513504E-02,.3493177E-02,.3472754E-02,
     $.3452237E-02,.3431624E-02,.3410918E-02,.3390119E-02,.3369227E-02,
     $.3348242E-02,.3327166E-02,.3305999E-02,.3284741E-02,.3263394E-02,
     $.3241957E-02,.3220431E-02,.3198818E-02,.3177116E-02,.3155328E-02,
     $.3133454E-02,.3111493E-02,.3089448E-02,.3067318E-02,.3045104E-02,
     $.3022806E-02,.3000426E-02,.2977964E-02,.2955420E-02,.2932796E-02,
     $.2910091E-02,.2887306E-02,.2864443E-02,.2841501E-02,.2818481E-02,
     $.2795384E-02,.2772211E-02,.2748961E-02,.2725637E-02,.2702238E-02/
      data (wp(i),i= 881, 920)/
     $.2678765E-02,.2655218E-02,.2631599E-02,.2607908E-02,.2584146E-02,
     $.2560312E-02,.2536409E-02,.2512437E-02,.2488395E-02,.2464286E-02,
     $.2440109E-02,.2415865E-02,.2391555E-02,.2367179E-02,.2342739E-02,
     $.2318235E-02,.2293667E-02,.2269037E-02,.2244344E-02,.2219590E-02,
     $.2194775E-02,.2169901E-02,.2144966E-02,.2119973E-02,.2094922E-02,
     $.2069814E-02,.2044649E-02,.2019428E-02,.1994152E-02,.1968821E-02,
     $.1943437E-02,.1917999E-02,.1892508E-02,.1866966E-02,.1841373E-02,
     $.1815729E-02,.1790036E-02,.1764294E-02,.1738503E-02,.1712665E-02/
      data (wp(i),i= 921, 960)/
     $.1686780E-02,.1660848E-02,.1634872E-02,.1608850E-02,.1582785E-02,
     $.1556676E-02,.1530525E-02,.1504331E-02,.1478097E-02,.1451822E-02,
     $.1425507E-02,.1399154E-02,.1372762E-02,.1346332E-02,.1319866E-02,
     $.1293363E-02,.1266825E-02,.1240253E-02,.1213646E-02,.1187006E-02,
     $.1160334E-02,.1133630E-02,.1106895E-02,.1080130E-02,.1053335E-02,
     $.1026511E-02,.9996593E-03,.9727801E-03,.9458743E-03,.9189426E-03,
     $.8919858E-03,.8650045E-03,.8379996E-03,.8109717E-03,.7839217E-03,
     $.7568502E-03,.7297579E-03,.7026457E-03,.6755143E-03,.6483644E-03/
      data (wp(i),i= 961,1000)/
     $.6211967E-03,.5940120E-03,.5668111E-03,.5395947E-03,.5123635E-03,
     $.4851182E-03,.4578597E-03,.4305887E-03,.4033058E-03,.3760120E-03,
     $.3487078E-03,.3213941E-03,.2940716E-03,.2667411E-03,.2394033E-03,
     $.2120589E-03,.1847087E-03,.1573535E-03,.1299941E-03,.1026314E-03,
     $.7526651E-04,.4790311E-04,.2057885E-04,.3924530E-02,.3924469E-02,
     $.3924348E-02,.3924167E-02,.3923925E-02,.3923623E-02,.3923260E-02,
     $.3922837E-02,.3922354E-02,.3921810E-02,.3921206E-02,.3920541E-02,
     $.3919816E-02,.3919030E-02,.3918185E-02,.3917278E-02,.3916312E-02/
      data (wp(i),i=1001,1040)/
     $.3915285E-02,.3914198E-02,.3913051E-02,.3911843E-02,.3910575E-02,
     $.3909247E-02,.3907858E-02,.3906410E-02,.3904901E-02,.3903332E-02,
     $.3901703E-02,.3900014E-02,.3898265E-02,.3896456E-02,.3894587E-02,
     $.3892658E-02,.3890668E-02,.3888619E-02,.3886510E-02,.3884342E-02,
     $.3882113E-02,.3879825E-02,.3877476E-02,.3875068E-02,.3872601E-02,
     $.3870074E-02,.3867487E-02,.3864840E-02,.3862134E-02,.3859369E-02,
     $.3856544E-02,.3853660E-02,.3850716E-02,.3847713E-02,.3844651E-02,
     $.3841530E-02,.3838349E-02,.3835109E-02,.3831811E-02,.3828453E-02/
      data (wp(i),i=1041,1080)/
     $.3825036E-02,.3821561E-02,.3818026E-02,.3814433E-02,.3810781E-02,
     $.3807070E-02,.3803300E-02,.3799473E-02,.3795586E-02,.3791641E-02,
     $.3787638E-02,.3783576E-02,.3779456E-02,.3775278E-02,.3771042E-02,
     $.3766747E-02,.3762395E-02,.3757984E-02,.3753516E-02,.3748990E-02,
     $.3744406E-02,.3739765E-02,.3735066E-02,.3730309E-02,.3725495E-02,
     $.3720624E-02,.3715695E-02,.3710709E-02,.3705666E-02,.3700566E-02,
     $.3695408E-02,.3690194E-02,.3684923E-02,.3679596E-02,.3674211E-02,
     $.3668770E-02,.3663273E-02,.3657719E-02,.3652109E-02,.3646442E-02/
      data (wp(i),i=1081,1120)/
     $.3640720E-02,.3634941E-02,.3629106E-02,.3623216E-02,.3617269E-02,
     $.3611267E-02,.3605209E-02,.3599096E-02,.3592927E-02,.3586703E-02,
     $.3580424E-02,.3574090E-02,.3567700E-02,.3561256E-02,.3554757E-02,
     $.3548203E-02,.3541594E-02,.3534931E-02,.3528213E-02,.3521441E-02,
     $.3514615E-02,.3507734E-02,.3500800E-02,.3493812E-02,.3486770E-02,
     $.3479674E-02,.3472524E-02,.3465321E-02,.3458065E-02,.3450756E-02,
     $.3443393E-02,.3435977E-02,.3428508E-02,.3420987E-02,.3413413E-02,
     $.3405786E-02,.3398107E-02,.3390375E-02,.3382592E-02,.3374756E-02/
      data (wp(i),i=1121,1160)/
     $.3366868E-02,.3358928E-02,.3350937E-02,.3342894E-02,.3334800E-02,
     $.3326654E-02,.3318457E-02,.3310208E-02,.3301909E-02,.3293559E-02,
     $.3285158E-02,.3276707E-02,.3268205E-02,.3259653E-02,.3251051E-02,
     $.3242398E-02,.3233696E-02,.3224944E-02,.3216142E-02,.3207290E-02,
     $.3198390E-02,.3189440E-02,.3180440E-02,.3171392E-02,.3162295E-02,
     $.3153149E-02,.3143955E-02,.3134712E-02,.3125421E-02,.3116082E-02,
     $.3106695E-02,.3097260E-02,.3087778E-02,.3078247E-02,.3068670E-02,
     $.3059045E-02,.3049373E-02,.3039654E-02,.3029888E-02,.3020075E-02/
      data (wp(i),i=1161,1200)/
     $.3010217E-02,.3000311E-02,.2990360E-02,.2980362E-02,.2970318E-02,
     $.2960229E-02,.2950094E-02,.2939914E-02,.2929688E-02,.2919418E-02,
     $.2909102E-02,.2898742E-02,.2888336E-02,.2877887E-02,.2867393E-02,
     $.2856855E-02,.2846273E-02,.2835647E-02,.2824977E-02,.2814264E-02,
     $.2803508E-02,.2792708E-02,.2781865E-02,.2770980E-02,.2760052E-02,
     $.2749081E-02,.2738068E-02,.2727013E-02,.2715915E-02,.2704776E-02,
     $.2693596E-02,.2682374E-02,.2671110E-02,.2659805E-02,.2648460E-02,
     $.2637073E-02,.2625646E-02,.2614179E-02,.2602671E-02,.2591123E-02/
      data (wp(i),i=1201,1240)/
     $.2579536E-02,.2567908E-02,.2556241E-02,.2544535E-02,.2532789E-02,
     $.2521005E-02,.2509181E-02,.2497319E-02,.2485419E-02,.2473480E-02,
     $.2461503E-02,.2449488E-02,.2437436E-02,.2425346E-02,.2413218E-02,
     $.2401054E-02,.2388852E-02,.2376614E-02,.2364339E-02,.2352028E-02,
     $.2339680E-02,.2327296E-02,.2314877E-02,.2302422E-02,.2289931E-02,
     $.2277405E-02,.2264844E-02,.2252249E-02,.2239618E-02,.2226953E-02,
     $.2214254E-02,.2201520E-02,.2188753E-02,.2175952E-02,.2163117E-02,
     $.2150249E-02,.2137349E-02,.2124415E-02,.2111448E-02,.2098449E-02/
      data (wp(i),i=1241,1280)/
     $.2085417E-02,.2072354E-02,.2059258E-02,.2046131E-02,.2032972E-02,
     $.2019782E-02,.2006561E-02,.1993309E-02,.1980026E-02,.1966713E-02,
     $.1953370E-02,.1939996E-02,.1926592E-02,.1913159E-02,.1899697E-02,
     $.1886205E-02,.1872684E-02,.1859134E-02,.1845555E-02,.1831949E-02,
     $.1818314E-02,.1804650E-02,.1790960E-02,.1777241E-02,.1763495E-02,
     $.1749722E-02,.1735922E-02,.1722096E-02,.1708242E-02,.1694363E-02,
     $.1680457E-02,.1666526E-02,.1652569E-02,.1638586E-02,.1624578E-02,
     $.1610545E-02,.1596488E-02,.1582405E-02,.1568299E-02,.1554168E-02/
      data (wp(i),i=1281,1320)/
     $.1540013E-02,.1525835E-02,.1511633E-02,.1497407E-02,.1483159E-02,
     $.1468888E-02,.1454594E-02,.1440278E-02,.1425940E-02,.1411579E-02,
     $.1397197E-02,.1382794E-02,.1368369E-02,.1353923E-02,.1339456E-02,
     $.1324969E-02,.1310461E-02,.1295933E-02,.1281385E-02,.1266817E-02,
     $.1252230E-02,.1237623E-02,.1222998E-02,.1208353E-02,.1193690E-02,
     $.1179009E-02,.1164309E-02,.1149592E-02,.1134857E-02,.1120104E-02,
     $.1105334E-02,.1090547E-02,.1075743E-02,.1060923E-02,.1046086E-02,
     $.1031234E-02,.1016365E-02,.1001481E-02,.9865808E-03,.9716658E-03/
      data (wp(i),i=1321,1360)/
     $.9567359E-03,.9417913E-03,.9268321E-03,.9118587E-03,.8968712E-03,
     $.8818700E-03,.8668551E-03,.8518269E-03,.8367855E-03,.8217313E-03,
     $.8066644E-03,.7915851E-03,.7764936E-03,.7613902E-03,.7462750E-03,
     $.7311483E-03,.7160104E-03,.7008614E-03,.6857017E-03,.6705313E-03,
     $.6553507E-03,.6401600E-03,.6249593E-03,.6097491E-03,.5945295E-03,
     $.5793007E-03,.5640630E-03,.5488166E-03,.5335618E-03,.5182987E-03,
     $.5030277E-03,.4877489E-03,.4724626E-03,.4571690E-03,.4418684E-03,
     $.4265610E-03,.4112470E-03,.3959267E-03,.3806003E-03,.3652680E-03/
      data (wp(i),i=1361,1400)/
     $.3499301E-03,.3345867E-03,.3192383E-03,.3038849E-03,.2885269E-03,
     $.2731644E-03,.2577977E-03,.2424270E-03,.2270526E-03,.2116747E-03,
     $.1962935E-03,.1809093E-03,.1655224E-03,.1501328E-03,.1347410E-03,
     $.1193471E-03,.1039514E-03,.8855408E-04,.7315545E-04,.5775582E-04,
     $.4235569E-04,.2695689E-04,.1158044E-04,.3140018E-02,.3139987E-02,
     $.3139926E-02,.3139833E-02,.3139709E-02,.3139554E-02,.3139368E-02,
     $.3139152E-02,.3138904E-02,.3138625E-02,.3138316E-02,.3137975E-02,
     $.3137604E-02,.3137201E-02,.3136768E-02,.3136304E-02,.3135809E-02/
      data (wp(i),i=1401,1440)/
     $.3135283E-02,.3134726E-02,.3134138E-02,.3133519E-02,.3132869E-02,
     $.3132189E-02,.3131477E-02,.3130735E-02,.3129962E-02,.3129158E-02,
     $.3128323E-02,.3127457E-02,.3126560E-02,.3125633E-02,.3124675E-02,
     $.3123686E-02,.3122666E-02,.3121615E-02,.3120534E-02,.3119422E-02,
     $.3118279E-02,.3117105E-02,.3115901E-02,.3114666E-02,.3113400E-02,
     $.3112103E-02,.3110776E-02,.3109418E-02,.3108029E-02,.3106610E-02,
     $.3105160E-02,.3103680E-02,.3102169E-02,.3100627E-02,.3099055E-02,
     $.3097452E-02,.3095819E-02,.3094155E-02,.3092461E-02,.3090736E-02/
      data (wp(i),i=1441,1480)/
     $.3088981E-02,.3087195E-02,.3085379E-02,.3083532E-02,.3081655E-02,
     $.3079748E-02,.3077810E-02,.3075842E-02,.3073843E-02,.3071815E-02,
     $.3069756E-02,.3067666E-02,.3065547E-02,.3063397E-02,.3061217E-02,
     $.3059007E-02,.3056766E-02,.3054496E-02,.3052195E-02,.3049865E-02,
     $.3047504E-02,.3045113E-02,.3042692E-02,.3040242E-02,.3037761E-02,
     $.3035250E-02,.3032709E-02,.3030139E-02,.3027538E-02,.3024908E-02,
     $.3022248E-02,.3019558E-02,.3016838E-02,.3014089E-02,.3011310E-02,
     $.3008501E-02,.3005662E-02,.3002794E-02,.2999896E-02,.2996969E-02/
      data (wp(i),i=1481,1520)/
     $.2994012E-02,.2991026E-02,.2988010E-02,.2984965E-02,.2981890E-02,
     $.2978786E-02,.2975652E-02,.2972489E-02,.2969297E-02,.2966075E-02,
     $.2962825E-02,.2959545E-02,.2956236E-02,.2952897E-02,.2949530E-02,
     $.2946134E-02,.2942708E-02,.2939254E-02,.2935770E-02,.2932258E-02,
     $.2928716E-02,.2925146E-02,.2921547E-02,.2917919E-02,.2914262E-02,
     $.2910577E-02,.2906863E-02,.2903120E-02,.2899349E-02,.2895549E-02,
     $.2891720E-02,.2887863E-02,.2883978E-02,.2880064E-02,.2876122E-02,
     $.2872151E-02,.2868152E-02,.2864125E-02,.2860069E-02,.2855985E-02/
      data (wp(i),i=1521,1560)/
     $.2851873E-02,.2847734E-02,.2843565E-02,.2839369E-02,.2835145E-02,
     $.2830893E-02,.2826613E-02,.2822305E-02,.2817970E-02,.2813606E-02,
     $.2809215E-02,.2804796E-02,.2800350E-02,.2795875E-02,.2791374E-02,
     $.2786844E-02,.2782288E-02,.2777704E-02,.2773092E-02,.2768453E-02,
     $.2763787E-02,.2759094E-02,.2754373E-02,.2749625E-02,.2744850E-02,
     $.2740048E-02,.2735219E-02,.2730363E-02,.2725480E-02,.2720571E-02,
     $.2715634E-02,.2710671E-02,.2705681E-02,.2700664E-02,.2695621E-02,
     $.2690551E-02,.2685454E-02,.2680331E-02,.2675182E-02,.2670006E-02/
      data (wp(i),i=1561,1600)/
     $.2664804E-02,.2659576E-02,.2654321E-02,.2649040E-02,.2643733E-02,
     $.2638400E-02,.2633041E-02,.2627657E-02,.2622246E-02,.2616809E-02,
     $.2611347E-02,.2605858E-02,.2600344E-02,.2594805E-02,.2589240E-02,
     $.2583649E-02,.2578033E-02,.2572391E-02,.2566724E-02,.2561032E-02,
     $.2555315E-02,.2549572E-02,.2543804E-02,.2538011E-02,.2532193E-02,
     $.2526350E-02,.2520483E-02,.2514590E-02,.2508672E-02,.2502730E-02,
     $.2496763E-02,.2490772E-02,.2484756E-02,.2478715E-02,.2472650E-02,
     $.2466561E-02,.2460447E-02,.2454309E-02,.2448147E-02,.2441961E-02/
      data (wp(i),i=1601,1640)/
     $.2435751E-02,.2429516E-02,.2423258E-02,.2416976E-02,.2410670E-02,
     $.2404340E-02,.2397986E-02,.2391609E-02,.2385209E-02,.2378784E-02,
     $.2372337E-02,.2365866E-02,.2359371E-02,.2352853E-02,.2346312E-02,
     $.2339748E-02,.2333161E-02,.2326551E-02,.2319918E-02,.2313262E-02,
     $.2306584E-02,.2299882E-02,.2293158E-02,.2286411E-02,.2279642E-02,
     $.2272850E-02,.2266036E-02,.2259200E-02,.2252341E-02,.2245460E-02,
     $.2238557E-02,.2231631E-02,.2224684E-02,.2217715E-02,.2210724E-02,
     $.2203711E-02,.2196677E-02,.2189620E-02,.2182543E-02,.2175443E-02/
      data (wp(i),i=1641,1680)/
     $.2168323E-02,.2161180E-02,.2154017E-02,.2146832E-02,.2139626E-02,
     $.2132400E-02,.2125152E-02,.2117883E-02,.2110593E-02,.2103282E-02,
     $.2095951E-02,.2088599E-02,.2081226E-02,.2073833E-02,.2066420E-02,
     $.2058986E-02,.2051531E-02,.2044057E-02,.2036562E-02,.2029047E-02,
     $.2021513E-02,.2013958E-02,.2006384E-02,.1998789E-02,.1991175E-02,
     $.1983542E-02,.1975888E-02,.1968216E-02,.1960524E-02,.1952812E-02,
     $.1945082E-02,.1937332E-02,.1929563E-02,.1921775E-02,.1913968E-02,
     $.1906142E-02,.1898298E-02,.1890434E-02,.1882552E-02,.1874652E-02/
      data (wp(i),i=1681,1720)/
     $.1866733E-02,.1858795E-02,.1850840E-02,.1842866E-02,.1834874E-02,
     $.1826863E-02,.1818835E-02,.1810789E-02,.1802725E-02,.1794643E-02,
     $.1786544E-02,.1778427E-02,.1770292E-02,.1762140E-02,.1753970E-02,
     $.1745784E-02,.1737580E-02,.1729359E-02,.1721120E-02,.1712865E-02,
     $.1704593E-02,.1696304E-02,.1687999E-02,.1679677E-02,.1671338E-02,
     $.1662983E-02,.1654611E-02,.1646223E-02,.1637819E-02,.1629399E-02,
     $.1620962E-02,.1612510E-02,.1604042E-02,.1595557E-02,.1587058E-02,
     $.1578542E-02,.1570011E-02,.1561465E-02,.1552903E-02,.1544325E-02/
      data (wp(i),i=1721,1760)/
     $.1535733E-02,.1527125E-02,.1518503E-02,.1509865E-02,.1501213E-02,
     $.1492545E-02,.1483863E-02,.1475167E-02,.1466456E-02,.1457730E-02,
     $.1448990E-02,.1440236E-02,.1431467E-02,.1422684E-02,.1413888E-02,
     $.1405077E-02,.1396253E-02,.1387414E-02,.1378563E-02,.1369697E-02,
     $.1360818E-02,.1351926E-02,.1343020E-02,.1334101E-02,.1325169E-02,
     $.1316224E-02,.1307265E-02,.1298294E-02,.1289310E-02,.1280314E-02,
     $.1271305E-02,.1262283E-02,.1253249E-02,.1244202E-02,.1235143E-02,
     $.1226072E-02,.1216989E-02,.1207894E-02,.1198787E-02,.1189668E-02/
      data (wp(i),i=1761,1800)/
     $.1180538E-02,.1171396E-02,.1162242E-02,.1153077E-02,.1143900E-02,
     $.1134712E-02,.1125513E-02,.1116303E-02,.1107082E-02,.1097850E-02,
     $.1088607E-02,.1079354E-02,.1070089E-02,.1060815E-02,.1051529E-02,
     $.1042234E-02,.1032928E-02,.1023612E-02,.1014286E-02,.1004950E-02,
     $.9956034E-03,.9862475E-03,.9768819E-03,.9675067E-03,.9581219E-03,
     $.9487276E-03,.9393240E-03,.9299112E-03,.9204892E-03,.9110581E-03,
     $.9016180E-03,.8921690E-03,.8827112E-03,.8732448E-03,.8637697E-03,
     $.8542861E-03,.8447941E-03,.8352937E-03,.8257851E-03,.8162684E-03/
      data (wp(i),i=1801,1840)/
     $.8067436E-03,.7972109E-03,.7876703E-03,.7781220E-03,.7685659E-03,
     $.7590023E-03,.7494312E-03,.7398528E-03,.7302670E-03,.7206740E-03,
     $.7110739E-03,.7014668E-03,.6918528E-03,.6822320E-03,.6726045E-03,
     $.6629703E-03,.6533295E-03,.6436824E-03,.6340289E-03,.6243691E-03,
     $.6147032E-03,.6050312E-03,.5953533E-03,.5856694E-03,.5759799E-03,
     $.5662846E-03,.5565837E-03,.5468774E-03,.5371657E-03,.5274486E-03,
     $.5177264E-03,.5079991E-03,.4982667E-03,.4885295E-03,.4787874E-03,
     $.4690406E-03,.4592892E-03,.4495332E-03,.4397729E-03,.4300082E-03/
      data (wp(i),i=1841,1880)/
     $.4202392E-03,.4104661E-03,.4006890E-03,.3909079E-03,.3811229E-03,
     $.3713342E-03,.3615418E-03,.3517459E-03,.3419465E-03,.3321437E-03,
     $.3223377E-03,.3125285E-03,.3027162E-03,.2929009E-03,.2830827E-03,
     $.2732617E-03,.2634380E-03,.2536118E-03,.2437830E-03,.2339519E-03,
     $.2241184E-03,.2142827E-03,.2044449E-03,.1946051E-03,.1847634E-03,
     $.1749198E-03,.1650745E-03,.1552276E-03,.1453792E-03,.1355293E-03,
     $.1256781E-03,.1158257E-03,.1059721E-03,.9611747E-04,.8626190E-04,
     $.7640548E-04,.6654832E-04,.5669051E-04,.4683217E-04,.3697344E-04/
      data (wp(i),i=1881,1883)/
     $.2711461E-04,.1725677E-04,.7413338E-05/
@}

\begin{thebibliography}{99}
\bibitem{GeorgeandLiu}
George, A. and Liu, J. W. H;
``Computer Solution of Sparse Positive Definite Systems'',
Prentice-Hall Inc, 1981

\end{thebibliography}
\end{document}
