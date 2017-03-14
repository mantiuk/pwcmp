# pwcmp 

This is a set of matlab functions for scaling of pairwise comparison experiment results based on the Thurstone's model V assumptions.

The main features:

* The scaling can work with imbalanced and incomplete data, in which not all pairs are compared and some pairs are compared more often than the others.

* Additional priors reduce bias due to the non-linear nature of the problem.

* The code can compute confidence intervals using bootstapping. 

## Usage

See examples in "examples" directory.

## Revision history

v0.1 (14.12.2016) Internal release

## Literature

There is a number of papers describing the technique. A similar technique has been described in:

D. Silverstein and J. Farrell, “Efficient method for paired comparison,” J. Electron. Imaging, vol. 10, no. 2, pp. 394–398, 2001.

However, the code contains a few improvements that make scaling more
robust and less prone to bias. Those extensions will be documented in
a separate report (to be published).

## License

MIT
