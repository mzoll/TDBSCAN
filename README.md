TDBSCAN: Temporal-Density-Based Spatial Clustering for Applications with Noise
===

# Introduction
TDBScan is an portable algorithm that allows the identification of clusters of causally connected blibs, commonly referred 
to as Signal, in a sea of other quasi-randomly occurring blibs, commonly referred to as noise, which do not share a common
generation process, aka. are not from the same source.

The algorithm works on a time-series of blibs, where a _Blib_ is one occurrence, event or information point in time.
In this context, a __Blib__ thus consists of a time-index and an Ordinate. The ordinate part may be any sort of data or 
set of variables; in a simple case it might be just be represented by a single continuous variable, which can be seen as
a location. We will later see that the Ordinate can also take the shape of a 3d-vector in space, thus relating the scenario
to applications in the real world. It might, however, also be something more complex with a variety of varibales:
contiuous, descrete or even categorical. The better one understands the relations and mechanisms governing these varibales and also their evolment and entwinement with time, 
the better one will be able to the set approprate parameters for the algorithm. In last consequence the better one is able to describe and model 
the mechnism driving the creation of said blibs as being caused by the Signal source, one can really make the algorithm excell.

The algorithm can be driven with a fixed set of blibs fed to it, or can be run in a continious mode feeding it new blibs bite for bite 
as they become available in real-time.

The algorithm shares many ideas and traits with its name-giving brother algorithm __DBScan__. However, by the utilization 
of the monotonic continuous variable time, it is able to run very effectively through vast data-sets and keep memory footprint and cpu utilization at a minimum. 

# Applications
As a theoretical problem setting, and very similar to its namesake, the algorithm can be allied to cluster blibs in fixed
size datasets in applications of data sciences, which are either artificially constructed or already are an abstraction themselves.
But the algorithm can also run on real-life problems: 
One of the most prominent are segemented scientific instruments, __detectors__, which continioulsy monitor
for occuences of scientific interest. If such an incident should occure they record the morphology of such an event and 
from there derive properties that allow inference of the physical cause.
Such detectors most often consist of many isolated subdetectors instrumenting a single probe volume, having a limit set of states and readouts:
either simple on/off, thus detection/non-detection, or they have a quantive dimension to the detection (signal-level, voltage, number of photons/elektrons, ...).
However, in real-world applications the subdetectors always suffer from sources of noisance disturbance, so called noise,
causing a __triggering__ of the detection-readout where no signal source had been present in the probe volume.
Noise can be the dominant part of subdetector-triggers in some experiments, thus effective ways to identify those
occurences where an acutal signal source was present in the detector assembly, and rejecting all other noise triggers is the first step.

TDBScan can do just that. By bottom-up modelling the properties of the physical signal source, its 
morphology and time-development, as well as taking into account the detector setup, for example its spatial distribution,
and detector efficiency and response, noise triggers can be very effectively separated from those actually caused by signal.

In application of the algorithm it is not important to exactly physically model the complete causal chain by means of
translation and propagation in time. It suffices to provide an idea of how detector triggers that are close in time relate to each other.

## IceCube
This algorithm was developed with one such detector in mind, the neutrino observatory IceCube, which is build at the South Pole.
It consists of many thousands of photon detectors (Photomultiplier-tubes) frozen in ice in a regular pattern, instrumenting a cubic kilometre of probe volume. Neutrinos
traversing the detector and undergoing interaction, will be seen as flashes of light (typical length of micro seconds) of specific morphology within the probe volume.
In the occurrence of such event, detectors in the vicinity of such event can be hit by those photons, trigger, and send their readout to the data collection on surface.
The simplest abstraction of such read-out is a __Hit__, which consist merely of a time-stamp and the location of the subdetector. 
The pattern in time as these hits are developing in the probe volume allow the interferrence of interaction type the neutrino underwent in the first place and its energy and direction.
However, the technology of the photomuliplier tubes is suffering from thermal emission of electrons along the same detection chain,
leading to said nuisance triggers, which is the primary source of noise in the detector. Among other processes in the detector, 
which do not relate to events caused neutrinos, the level of noisance or noise hits in the detector outweights those caused by actual signal by a factor 1000 <span style="color: red;"> (from the top of my head, needs checking!!!) </span>.

* (MORE!)
* The algorithm was applied ...
* Successfully run on raw detector readout from the IceCube Detector in 86-string configuration.
* Raw-detector triggers frequency ~Ghz (SUM of all DOMs)
* Running at near-realtime on a single CPU
* Link Licenciate thesis
* Link PhD thesis


# Working principle

## Abstract Description
This is a very top level description about the algorithm:

At its core the TDBSCAN Algorithm works as a time progressive sorting machine, attributing blibs to one or more already identified 
Clusters of blibs:

By feeding the __next blib__ ($b_i(t_{i},)$) in time (time-order!) into the machine, these clusters ($K$) will form, evolve, but
also eventually die; which is when no more blibs can be added to them. Thus at all times, there is a varying number of Clusters, 
which differ in shape and size, what we can refer to as morphology. At a minimum there is always one cluster present, which is the trivial Cluster,
formed by the previous blib ($b_{i-1}(t_{i-1}); t_{i-1} <= t_{i}$) by itself alone.

This current blib $b_{\text{c}}(t_\text{c},o)$, defining by its time coordinate the 'now' for the algorithm ($t_{\text{now}} = b_{\text{c}}(t_\text{c})$),
is upon it's processing evaluated by a user-defined but fixed
comparison prescription, the so called __Connector__ ($C(b_{n}, b_{c}); b_{n} \in K$), if it causally relates to the other blib within the cluster. 
If this is to be found true this forms __evidence__ that the blib is causally connected with the cluster. The algorithm requires that there is enough 
evidence collected before a blib is to be considered causally connected to the cluster as a whole. Currently, this is implemented
as a multiplicity requirement, so that at least $m_{\text{evid}}$
causal connections with blibs already in the cluster must be found. Only those blibs, which lie reasonably close in the past, 
described by a permitted time-window $t_{\text{evid}}$, are considered for this criterion.  
If the evidence reaches the required limit the blib is effectively added to the cluster.

As clusters in their early development (seeding), the so called __emerging__ clusters, with less blibs then $m_{\text{evid}}-1$ cannot fulfil this criterion,
a different but somewhat stronger criterion is applied requiring that all blibs are causally connected to the new blib. 
As more blibs are added, an previously emerging cluster eventually becomes an __active__ cluster once its multiplicity 
reaches the nominally required $m_{\text{evid}}$ blibs.

The time-window requirement results in that only the most recent blibs participate in the evaluation if a new blib should be added. 
In an abstract view this means that all blibs which are within this time-window form, the so called __active volume__ in the complex space formed
by varying variable time and ordinate. The morphology of this active volume exerted by its constituents and how the next blib relates to it, 
forms the driving principle for the further growth of any cluster.
Once having added a blib to a cluster it itself will be part of that active volume and influence its morphology when checking the next upcoming future blib.

With progression of time the addition of blibs to a cluster, its growth, can eventually stagnate. With more and more blibs falling out of the active volume as 
time moves forward, because they are then too far removed in time ($ b_{\text{past}}(t); b \in K < t_\text{now} - t_{\text{evid}}$), the active volume shrinks until
a point where less than $m_{\text{evid}}$ of the cluster's constituents are actively participating, and the evidence requirement for adding further blibs can no longer be fulfilled.
The cluster is then effectively concluded.
Some __emerging__ clusters might not even make it out of infancy before die off. Only clusters which once where __active__ will be presented as factual clusters. 

In extension that means that blibs, which never were added to any cluster making it to this adult stage, are regarded as 
not belonging to any factual cluster, and so are considered as noise.

In processing the algorithm keeps track of the blibs and clusters percolating through it, promoting clusters on the fly 
and pushing factual clusters to the output as soon as they deterministically concluded.
The output then is a list of sets of time-ordered blibs, each one being one factual cluster.

## Deep dive

