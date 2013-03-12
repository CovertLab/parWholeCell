parWholeCell
============
Project to develop a parallelized whole-cell model simulator. 

Currently implements a simplified whole-cell model using object-oriented MATLAB. 
![Model](https://github.com/CovertLab/parWholeCell/raw/master/doc/model.png)
![Simulation architecture](https://github.com/CovertLab/parWholeCell/raw/master/doc/architecture.png)

As of 3/12/2013 the code implements a simple whole-cell model containing 7 submodels. The submodels are extremely simple.
* Complexation
* Metabolism: flux-balance analysis
* RNA degradation
* RNA maturation: processing, cleavage, etc.
* Protein maturation: maturation, translocation
* Transcription
* Translation

and 3 states:
* MoleculeCounts: metabolite, RNA, protein monomer, complexes. Each RNA and protein has two forms: nascent and mature.
* Metabolism: growth and reaction fluxes
* Mass

The model includes 3 compartments:
* Cytosol
* Membrane
* Extracellular space


Requirements
-------------------------
* MATLAB >= R2012b
* MATLAB toolboxes
	* bioinformatics

Installation
-------------------------
1. Clone from github repository
	* https://github.com/CovertLab/parWholeCell.git
	* ssh://git@github.com:CovertLab/parWholeCell.git
2. Open MATLAB, change to parWholeCell directory
3. In addition, setup the MATLAB warnings and path at the beginning of each MATLAB session
	```matlab
	setWarnings();
	setPath();
	```
4. Thats it!


How to run simulation and plot
-------------------------
Use the code below to run a simulation with the default options and parameter values.

```matlab
%set warnings and path (need to execute once per session)
setWarnings();
setPath();

%simulate and plot
runSimulation(); %example use in wholecell.sim.SimulationTest.testRunSimulation
```


How to run tests
-------------------------
parWholeCell uses [MATLAB XUnit](http://www.mathworks.com/matlabcentral/fileexchange/22846-matlab-xunit-test-framework) to run unit tests. Use the code below to build fixtures and run tests.

```matlab
%set warnings and path (need to execute once per session)
setWarnings();
setPath();

%generate fixtures for tests so tests run faster
generateTestFixtures();

%Run tests. Tests are logged to out/test/results.xml in JUnit-style XML
runTests();
```


Simulation API
-------------------------
```matlab
import wholecell.sim.Simulation;

sim = Simulation()

%get options, parameters
options = sim.getOptions();
parameters = sim.getParameters();

%set options, parameters
options = struct('lengthSec', 10);
sim.setOptions(options);
sim.setParameters(parameters);

%get states, process
s = sim.getState(<stateId>)
p = sim.getProcess(<processId>)

%run simulation
loggers = {wholecell.sim.logger.Shell()};
sim.run(<optional cell array of loggers>);

%calculate initial conditions
sim.calcInitialConditions();

%calculate one time step
sim.evolveState()
```

Comparison to Karr et al., 2012 version
-------------------------
The main ideas of this model are identical to the 2012 version. The models differs in a few ways:
* The Metabolite, RNA, ProteinMonomer, and ProteinComplex states have been merged into one state called MoleculeCounts
* The RNA and protein maturation submodels have been merged into 2 submodels
* The model includes only three compartments. The terminal organelle and nucleoid compartments have been eliminated
* As of 3/12/2013 there are no submodels to represent DNA replication, cytokinesis, or protein decay
* The submodels are extremely simple
* There are only 2 RNA and protein forms: nascent and mature

Additionally, the implementation differs in a few ways:
* There is no database implementation of the knowledge base. Instead one class wholecell.kb.KnowledgeBase provides the same functionality as before. This class reads data from a Excel workbook.
* The names of the submodels, states, and methods, and properties have been modified in a few places for clarity.
* Most importantly, the process-state communication has been improved, keeping in mind whats needed to parallelize the implementation. We've introduced a new concept called a partition which represents the interface between a process and a state. 
	The partitions do several things. 
	* They encapsulate the mapping between states and proceses. 
	* They also isolate the submodels such that they don't read/write the same piece of memory. 
	
	States have a list of partitions and several methods for managing them.
	* addPartition: called during model construction to create a new process-state mapping
	* partition: splits a state in substates/partitions assigned to each process
	* merge: merges substates back into the full state
* It is now possible to vary the time step size
* There is no master table of parameters yet. This is easy to implement. I just haven't done it yet.
* There's no code for model fitting or analysis.
* The disk logger has been simplified. The logging functionality is the same. But there no reindexing and less support for efficient retrieval.	

Credits
-------------------------
parWholeCell was developed by researchers at Stanford University and Intel:
* [Markus Covert Systems Biology Lab](http://covertlab.stanford.edu/), [Stanford University](http://www.stanford.edu)
	* [Markus Covert](http://bioengineering.stanford.edu/faculty/covert.html). Contact [mcovert@stanford.edu](mailto:mcovert@stanford.edu).
	* [Jonathan Karr](http://www.stanford.edu/~jkarr/). Contact [jkarr@stanford.edu](mailto:jkarr@stanford.edu).	
* [Parallel Computing Lab](http://pcl.intel-research.net), [Intel Corporation](http://www.intel.com/)
	* [Sanchit Misra](http://pcl.intel-research.net/people/sanchit.htm). Contact [sanchit.misra@intel.com](mailto:sanchit.misra@intel.com).
	* [Kiran Pamnany](http://pcl.intel-research.net/people/kiran.htm). Contact [kiran.pamnany@intel.com](mailto:kiran.pamnany@intel.com).
	* [Mikhail Smelyanskiy](http://pcl.intel-research.net/people/misha.htm). Contact [mikhail.smelyanskiy@intel.com](mailto:mikhail.smelyanskiy@intel.com)
	