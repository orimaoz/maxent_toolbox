<h2>maxent_toolbox</h2>


<p>Maximum entropy toolbox for MATLAB is a free, open-source toolbox for finding the maximum entropy distribution of training data, based on a set of constraints or observables over the data.</p> 

<p>
Maximum entropy models give the mathematically minimal probabilistic models of the states or configurations of a systems, given the mean values of some set of observed functions. Since the entropy of a distribution measures the randomness or lack of interaction among different variables, the minimally structured distribution given a set of observables is the distribution with the maximal entropy that is consistent with these observables. These models are used approximate probability distributions of discrete states from limited data, where standard approaches (such as counting the states) cannot be applied.
</p>

<p>This toolbox is intended for learning probability distributions of <b>binary activity patterns</b>, i.e. patterns in the form
1000110100. Examples of problems where one might want to learn such distributions are characterising the joint
activity of neural populations (0/1 denoting if a certain neuron fired within a small time window) or gene expression profiles
(0/1 denoting if a certain gene is being expressed). The toolbox takes as an input a set of samples of activity patterns and learns a model of the probability over all states, thus extrapolating to the entire distribution over all possible activity patterns.</p>

<p>
Mathematically, in its discrete form, if x<sub>i</sub>
 are the elements of the system (here variables taking discrete values), then the maximum entropy model for p(x<sub>1</sub>,x<sub>2</sub>…x<sub>n</sub>)
 which is consistent with the average values of a set of functions f<sub>1</sub>(x)...f<sub>k</sub>(x) has a unique solution in the form:
</p>

<p>
<img src="http://latex.codecogs.com/svg.latex?\hat{p}(x)=\exp[\sum_{i=1}^{k}\lambda_i{f_{i}(\vec{x})}])" border="0"/>
</p>

<p>where &#955;<sub>1</sub>...&#955;<sub>k</sub> are the model parameters (or Lagrange multipliers) which can be found numerically.</p>

<p>The user can choose between several variants of maximum entropy models, each relying on a different set of observables or constraints. The maximum entropy models currently supported by this package are
<ul>
<li>Independent model, which uses &#60;x<sub>i</sub>&#62;<sub>data</sub> as constraints</li>
<li><a href="http://www.weizmann.ac.il/neurobiology/labs/schneidman/Publications/Schneidman+al_2006-Nature.pdf">Pairwise maximum entropy model</a>: &#60;x<sub>i</sub>&#62;<sub>data</sub> and &#60;x<sub>i</sub>x<sub>j</sub>&#62;<sub>data</sub></li>
<li><a href="http://wp.ist.ac.at/group_tkacik/wp-content/uploads/2011/02/JStatMech_Simplest_Maxent.pdf">K-synchrony model</a>: &#60;&Sigma;<sub>i</sub>x<sub>i</sub>&#62;<sub>data</sub></li>
<li><a href="http://wp.ist.ac.at/group_tkacik/wp-content/uploads/2011/02/PLOSCompBio_KPairwise.pdf">K-Pairwise model</a></li>
<li>Arbitrary set of high-order correlations</li>
<li><a href="https://www.biorxiv.org/node/143383.abstract">RP (Random Projection) model</a></li>
<li>Any combination of the above models</li>
</ul>
</p>

<p>The toolbox automatically switches between exhaustive solutions for small (&#60;30) groups of variables and Markov Chain Monte Carlo (MCMC) methods for larger groups and can be used to learn distributions of up to several hundreds of binary variables. The software is provided as an installable toolbox for MATLAB, and most of the code is written in heavily optimized C++ precompiled for Windows (64 bit), OS X and Linux (CentOS).</p>
<h2>Download</h2>
<ul>
<li><a href="https://github.com/orimaoz/maxent_toolbox/releases/latest">Download latest version</a> packaged as an installable MATLAB toolbox with mex files precompiled for Windows 64bit, MacOS and Linux (CentOS 64bit).</li>
<li><a href="https://github.com/orimaoz/maxent_toolbox">Download most recent source code</a>.</li>
</ul>

<h2>Install</h2>
<p>Click on the downloaded .mltbx file to install it as a toolbox for MATLAB. To uninstall, open MATLAB and navigate to <i>Add-Ons -> Manage Add-Ons.</i></p>
        
<h2><a id="documentation" class="anchor" href="#documentation-list" aria-hidden="true"><span aria-hidden="true" class="octicon octicon-link"></span></a>Documentation</h2>
<p>If this is the first time you are using the toolbox, start by reading the <a href="quickstart.html">quick start guide</a> or take a look at some <a href="maxent_example.html">example code</a> (which you can also run by typing "maxent_example" after installing the toolbox).</p>
<p>Experienced users can refer to the <a href="function_reference.html">function reference</a> for details about each of the available functions in the toolbox.</p>
<p>You can also check out answers to some<a href="faq.html"> frequently asked questions</a>.</p>


<h2><a id="build" class="anchor" href="#build" aria-hidden="true"><span aria-hidden="true" class="octicon octicon-link"></span></a>Build</h2>
<p><a href="build.html">Build instructions</a> for Windows, MacOS and Unix. Note that the toolbox installer already contains precompiled binaries, so building the toolbox is only necessary if you use a nonstandard operating system or want to make changes in the code.</p>

<h2>Cite</h2>
<p>If you use this toolbox as part of a published academic work, please cite it as:</p>
<pre>
Ori Maoz and Elad Schneidman. maxent_toolbox: Maximum Entropy Toolbox for MATLAB, version 1.0.2. 2017. doi: 10.5281/zenodo.191625. URL https://orimaoz.github.io/maxent_toolbox.
</pre>
<p>The corresponding BibTeX citation:</p>
<pre>
@misc{Maoz2017_191625,
author = {Maoz, Ori and Schneidman, Elad},
title = {maxent{\_}toolbox: Maximum Entropy Toolbox for MATLAB, version 1.0.2},
year = {2017},
version = {1.02},
publisher = {Zenodo},
doi = {10.5281/zenodo.191625},
url = {https://orimaoz.github.io/maxent_toolbox}
}
</pre>


<h2>License</h2>
<p>The Maximum Entropy Toolbox for MATLAB is free and open source, distributed under the <a href="https://opensource.org/licenses/MIT">MIT license</a>.</p>

<h2>
<a id="authors-and-contributors" class="anchor" href="#authors-and-contributors" aria-hidden="true"><span aria-hidden="true" class="octicon octicon-link"></span></a>Authors and Contributors</h3>
<p>Ori Maoz, <a href="mailto:orimaoz@gmail.com"> orimaoz@gmail.com</a></p>
