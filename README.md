TDBSCAN: Temporal Density-Based Spatial Clustering for Applications with Noise
===
![image](./docs/resources/tdbscan_logo_white_1280x640.svg)

# Introduction
TDBScan is an portable algorithm that allows the identification of clusters of causally connected blibs, commonly referred 
to as Signal, in a sea of other quasi-randomly occurring blibs, commonly referred to as noise, which do not share a common
generation process, aka do not stem from the same source.

The algorithm works on a time-series of blibs, where a _Blib_ is one occurrence, event or point of information in time.
In this context, a _Blib_ simply consists of a time-index and an Ordinate. The ordinate part may be any sort of data or 
set of variables; in a simple case it might be just be represented by a single continuous variable, which can be seen as
a location. We will later see that the Ordinate can also take the shape of a 3d-vector in space, thus relating the scenario
to applications in the real world. It might, however, also be something more complex with a variety of varibales:
contiuous, descrete or even categorical. The better one understands the relations and mechanisms governing these varibales and also their evolment in time, 
the better one will be able to the set approprate parameters for the algorithm. In last consequence the better one is able to describe and model 
the mechnism driving the creation of said blibs as being caused by the _Signal_ source, one can really make the algorithm excell.

The algorithm can be applied with a fixed set of blibs fed to it, or can be run in a continuous mode feeding it new blibs bite for bite 
as they become available in real-time.

The algorithm shares many ideas and traits with its name-giving brother algorithm **DBScan**. However, by the exploitation 
of the monotonic continuous variable time, it is able to run very effectively through vast data-sets and keep memory footprint and cpu utilization at a minimum. 

# Applications
As a theoretical problem setting, and very similar to its namesake, the algorithm can be allied to cluster blibs in fixed
size datasets in applications of data sciences, which are either artificially constructed or already are an abstraction themselves.
But the algorithm can also run on real-life problems: 
One of the most prominent are segemented scientific instruments, **detectors**, which continioulsy monitor
for occuences of deeper scientific interest. If such an incident should occure they record the morphology of such an event and 
from there derive properties that allow inference of a physical cause.

Such detectors most often consist of many isolated subdetectors instrumenting a single probe volume, having a limit set of states and readouts:
either simple on/off, thus detection/non-detection, or they have a quantitative dimension to the detection (signal-level, voltage, number of photons/elektrons, ...).
However, in real-world applications the subdetectors always suffer from sources of nuisance disturbance, so-called _Noise_,
causing a **triggering** of the detection-readout where no signal source had been present in the probe volume.
In some scientific experiments noise can be the dominant part of subdetector-triggers. The first step is always to 
discriminate noise at the level of the subdetector itself. But the possibilities for this are limited, as most often
the detector triggers caused by (irreducible) noise and those caused by actual signal are inherently the same. 
Thus effective ways to further supress noise, and identify those incidents in time, where an actual signal source was
present in the whole of the detector assembly, and rejecting all other noise triggers, is the next step in the process.

TDBScan can do just that. By bottom-up modeling the properties of the physical signal source, its 
morphology and time-development, as well as taking into account the detector assembly's physical setup,
detection efficiency and response, noise triggers can be very effectively separated from those actually caused by signal sources.

In application of the algorithm it is not important to exactly physically model the complete causal chain by means of
translation and propagation in time. It suffices to provide just an idea of how subdetector triggers that are close in time relate to each other.

## IceCube
This algorithm was developed with one such detector in mind, the neutrino observatory IceCube, which is build at the South Pole.
It consists of many thousands of photon detectors (Photomultiplier-tubes) frozen in ice in a regular pattern, instrumenting a cubic kilometre of probe volume.

Neutrinos traversing the detector and undergoing interaction, will be seen as flashes of light (typical length of micro seconds) of specific morphology within the probe volume.
In the occurrence of such event, detection units in the vicinity of such event can be hit by those photons, trigger, and send their readout to the data collection on the ice's surface.
The simplest abstraction of such read-out is a **Hit**, which consist merely of a time-stamp and the location of the detection unit. 
The pattern in time as these hits are developing in the probe volume allow the interferrence of interaction type the neutrino underwent and its energy and direction.
However, the technology of the photomuliplier tubes is suffering from thermal emission of electrons along the same detection chain,
leading to said nuisance triggers, which is the primary source of noise in the detector. Among other processes in the detector, 
which do not relate to events caused neutrinos, the level of noisance or noise hits in the detector outweights those caused by 
actual signal by a factor 1000 <span style="color: red;"> (from the top of my head, needs checking!!!) </span>.

* (MORE!)
* The algorithm was applied ...
* Successfully run on raw detector readout from the IceCube Detector in 86-string configuration.
* Raw-detector triggers frequency ~Ghz (SUM of all DOMs)
* Running at near-realtime at filter level on single CPU (2015)
* Link Licenciate thesis: [Improved methods for solar Dark Matter searches with the IceCube neutrino telescope](http://www.diva-portal.org/smash/record.jsf?pid=diva2%3A891768) 
* Link PhD thesis: [A search for solar dark matter with the IceCube neutrino detector: Advances in data treatment and analysis technique](http://www.diva-portal.org/smash/record.jsf?pid=diva2%3A892438)
* The algorithm  used to derive result for the following paper: [Search for annihilating dark matter in the Sun with 3 years of IceCube data](https://arxiv.org/abs/1612.05949)
* Check out these videos about the algorithms performance on real-world data:
  * [Raw detector readout - no treatment](https://youtu.be/Kxfi-zDacxI)
  * [Raw detector readout, visually emph. clusters](https://youtu.be/8dxrlA95h_4)
  * [Clusters isolated, reconstruction applied](https://youtu.be/HGEuLnkxe3s)
    (NOTE: the particle-trajectory reconstruction 'line-objects' is not part of the here discussed algorithm)

## Performance now

The original algorithm employed with the IceCube project has emerged around 2011, but was then further developed and refined 
until 2015 into the project _IceHive_ for the specific application with IceCube data structures. 

Revisiting past work in August and September of 2026, the algorithmic portion of the _IceHive_ project was broken out, and rethought and finally code 
completely rewritten into what this project comprises. This allowed to make the algorithm more transparent in its working 
principle in code and generalize the framework for application in general data analysis of time-series.

In its current form, written in C++(17) the algorithm is: 

* Lightweight: No heavy dependencies to foreign libraries
* Fast (quantify this!)
* Efficient (minimal copying objects in memory, small memory footprint)
* Flexible (algorithm allows customization at parameter and functional level)
* Extensible (classes are all transparent)

# Working principle

## Abstract Description
This is a very top level description about the algorithm:

At its core the TDBSCAN Algorithm works as a time progressive sorting machine, attributing blibs to one or more already identified 
Clusters of blibs:

By feeding the **next blib** ($b_i(t_{i},)$) in time (time-order!) into the machine, these clusters ($K$) will form, evolve, but
also eventually die; which is when no more blibs can be added to them. Thus at all times, there is a varying number of Clusters, 
which differ in shape and size, what we can refer to as morphology. At a minimum there is always one cluster present, which is the trivial Cluster,
formed by the previous blib ($b_{i-1}(t_{i-1}); t_{i-1} <= t_{i}$) by itself alone.

This current blib $b_{\text{c}}(t_\text{c},o)$, defining by its time coordinate the 'now' for the algorithm ($t_{\text{now}} = b_{\text{c}}(t_\text{c})$),
is upon it's processing evaluated by a user-defined but fixed
comparison prescription, the so-called **Connector** ($C(b_{n}, b_{c}); b_{n} \in K$), if it causally relates to the other blib within the cluster. 
If this is to be found true this forms **evidence** that the blib is causally connected with the cluster. The algorithm requires that there is enough 
evidence collected before a blib is to be considered causally connected to the cluster as a whole. Currently, this is implemented
as a multiplicity requirement, so that at least $m_{\text{evid}}$
causal connections with blibs already in the cluster must be found. Only those blibs, which lie reasonably close in the past, 
described by a permitted time-window $t_{\text{evid}}$, are considered for this criterion.  
If the evidence reaches the required limit the blib is effectively added to the cluster.

As clusters in their early development (seeding), the so called **emerging** clusters, with less blibs then $m_{\text{evid}}-1$ cannot fulfil this criterion,
a different but somewhat stronger criterion is applied requiring that all blibs are causally connected to the new blib. 
As more blibs are added, an previously emerging cluster eventually becomes an **active** cluster once its multiplicity 
reaches the nominally required $m_{\text{evid}}$ blibs.

The time-window requirement results in that only the most recent blibs participate in the evaluation if a new blib should be added. 
In an abstract view this means that all blibs which are within this time-window form, the so called **active volume** in the complex space formed
by varying variable time and ordinate. The morphology of this active volume exerted by its constituents and how the next blib relates to it, 
forms the driving principle for the further growth of any cluster.
Once having added a blib to a cluster it itself will be part of that active volume and influence its morphology when checking the next upcoming future blib.

With progression of time the addition of blibs to a cluster, its growth, can eventually stagnate. With more and more blibs falling out of the active volume as 
time moves forward, because they are then too far removed in time ($ b_{\text{past}}(t); b \in K < t_\text{now} - t_{\text{evid}}$), the active volume shrinks until
a point where less than $m_{\text{evid}}$ of the cluster's constituents are actively participating, and the evidence requirement for adding further blibs can no longer be fulfilled.
The cluster is then effectively concluded.
Some **emerging** clusters might not even make it out of infancy before die off. Only clusters which once where **active** will be presented as factual clusters. 

In extension that means that blibs, which never were added to any cluster making it to this adult stage, are regarded as 
not belonging to any factual cluster, and so are considered as noise.

In processing the algorithm keeps track of the blibs and clusters percolating through it, promoting clusters on the fly 
and pushing factual clusters to the output as soon as they deterministically concluded.
The output then is a list of sets of time-ordered blibs, each one being one factual cluster.

## Deep dive

On the Implementation level the algorithm works like this:

The algorithm object holds a registry for each: emerging clusters, active clusters and concluded clusters.

0. The next fed Blib is inspected, its time stamp becomes the new 'now' time. Time-order of the fed Blibs is ensured.
1. Each preexisting **emerging cluster** is checked if all blibs contained within are still within the eligible evidence time window. 
   If not, the cluster is eliminated.
2. Each preexisting **emerging cluster** is checked if __all__ blibs can be causally connected via the Connector to the current Blib.
   If so the blib is added to that cluster. If the multiplicity of blibs within the cluster reaches the required multiplicity $m$, it is transferred to a 
   temporary __newly established cluster__ list.
3. The blib is by itself put into a __new__ cluster on the **emerging cluster** registry.
4. Each preexisting **active cluster** is checked if at least $m$ blibs remain within the eligible evidence time window. 
   If not the cluster is transferred to the registry for **concluded clusters**.
5. Each remaining **active cluster is checked if __evidence__ greater than the required multiplicity is found.
   If so the blib is added to the cluster.
6. Each cluster on the temporary __newly established cluster__ list is validated against all active clusters if their overlap in number of blibs is at
   least the multiplicity $m$ times a user adjustable ratio. This called the early merge mechanism. 
   If enough overlap is found, the __newly established cluster__ is absorbed into the preexisting active cluster.
7. Repeat from step 0. as long as there are more blibs; If no more blibs are expected, the algorithm is promted to **finalize**:
   * All remaining **emerging clusters** are eliminated.
   * All remaining **active clusters** are transferred to the registry for **concluded clusters**.
   * The registry of  concluded cluster is return as the output. 

### Implementation details

The code has been written, so that there is no fixed definition of a Blib, so that the user can customize it to his disgression.
Any Object can act as a Blib as long as it can provide the following traits:
* It needs a Time defined, that is retrievable: ``::getTime()``
* In comparison to another Blib, their difference in time needs to be defined: ``::timeTo(other)``
* It needs a Ordinate defined, that is retrievable: ``getOrdinate()``
* In comparison to another Blib, their distance given by their ordinates needs to be defined: ``::distanceTo(other)``
* Objects need to be id-comparable ``::operator==(other)``, as well as ordered ``::operator<(other)``  

Objects for Time and Ordinate themselves also need to be well defined in the mathematical sense, for example are well structured, ordered and so on.
There are abstract base classes provided which can be extended to the exact use-case at hand. However, most standard variable types (``int``, ``double``, ...),
already serve well as prototypes to construct these properties for simple cases of Time and Ordinate classes.
Blibs are provided through an abstract Proto-class as well.

To Provide some head-start the most common and generalized use-case classes are already defined and implemented for convenience:
* ``ScalarTime``: scalar (floating-point) time-like object with double precision, time differences also represented by a double.
* ``Position1d``: scalar (floating-point) ordinate-like object with double precision, differences also represented by a double.
* ``Position3d``: 3-tuple of x,y,z with double precision, the difference being the scalar-(dot)-product; in earnest the 3 vector known from school mathematics.

Also they combine in to the Blib classes
* ``ScalarBlib``: ``ScalarTime`` and ``Postion1d`` form a very simple formulation of a Blib
* ``Blib3d``: ``ScalarTime`` and ``Position3d`` form a very relatable real-world formulation of a Blib; this forms what 
  is oftimes called the 4-component space-time vector in physics.

If one does not want to implement custom Time or Ordinate objects the latter classes can also simply be extended by inheritance to arrive at a 
deeper (having further properties) but still algorithmic functional Blib object. 

## The Connector

(Description of how Connectors are to be envisioned and implemented)
