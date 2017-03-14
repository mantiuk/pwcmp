# pwcmp 

This is a set of matlab functions for scaling of pairwise comparison experiment results based on Thurstone's model V assumptions.

The main features:

* The scaling can work with imbalanced and incomplete data, in which not all pairs are compared and some pairs are compared more often than the others.

* Additional priors reduce bias due to the non-linear nature of the problem.

* The code can compute confidence intervals using bootstapping. 

## Usage

'''
% Simple example showing how to execute scaling method

% Comparison matrix. This is an incomplete design in which 1-st conditon was
% compared with the 2nd and 2nd with the 3rd. In both cases 75 observers
% seleted one condition and 25 the other.
D = [ 0   25  0;
      75  0  25;
      0   75 0  ];
   
Q = pw_scale( D );

display( Q )
'''

See more examples in "examples" directory.

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
