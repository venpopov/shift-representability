These documents describe an ongoing mathematical psychology project. I want you to design an implement a flexible simple and extensible system in R for numerically confirming the necessary and sufficient properties for a set of random variables or distributions to be considered shift-representable.

The theorems are stated algebraically, and for some important families they can be confirmed symbolically. However, for many a symbolic approach is not feasible. This is why we want to be to test sets of random variables with different distributions, and varying parameters of these distributions, whether computationally and numerically they satisfy the listed properties either exactly, or to get an estimate of how much they violate them, if they do.

There are two important classes of inputs to consider:

- distributions that can be evaluated with very high degree a precision that have standard reference implementations (such as pgamma, qgamma for the cdf and quantile function of the gamma distribution in R
- custom models that can be evaluated through simulation, and whose outputs are stochastic.

For the initial implementation you should focus on the easier first case, but still keep in mind the code structure to allow it to extend. 

- first read the attached files and determine what properties must be evaluated of distributions, and come up with a plan as to how to do it computationally.
- then consider different implementation options that woul allow flexible extensible work, as we might continuously add new distributions to test
- it should easy for the user to specify what distributions to include in a test, and to easily vary their parameters
- they should be able to do this for quick individual tests for iterative development
- they should also be able to batch test many distributions across a grid of parameters without having to do much additional scaffolding 
- when user specify distributions, they might have to provide the pdf / cdf / invcdf for each. or maybe use package "distributional" which can put distribution objects invectors and then directly apply cdf() or quantile() to the distribution objects
- think about the various analytical objects defined in the manuscript and what role they need to play in the solution

Your ultimate role is that of the architecture designer, and you should strive for simple but general solutions that will make the project much easier to work with and maintain
