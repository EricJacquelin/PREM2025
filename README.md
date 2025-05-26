The respository is related to a paper submitted to Probabilistic Engineering Mechanics in 2025. All the programs used to obtain the results presented in the paper can be found in this respository.

# 1\.Beginning
The main file that must be launched is 

\- **Either** Random\_Mode\_MAIN\_NI.m   
\- **Or**     Random\_Mode\_MAIN\_Intrusive.m  

depending on the intrusive/non-intrusive process you want to test.


# 2\.    Data
First the data are called in 

data\_NdofSystem.m
# 3\.  Direct MCS
\- Then a MCS is performed to get the random modes:

\- Random\_Mode\_MCS\_direct.m

Outputs: 

- Eigenvalues: Om2\_u(:,i\_u): squared eigenfrequency vector obtained for the i\_u-th set of uncertain parameters.
- Eigenvectors: Psi\_u(i\_mod).Psi\_u(:,i\_u): i\_mod-th eigenvector obtained for the i\_u-th set of uncertain parameters.
- M\_mod\_u (i\_mod,i\_u): modal mass vector obtained for the i\_u-th set of uncertain parameters: they are identified to build the FRF.


# 4\. Metamodel
Then a PCE metamodel is built. 

At that level, there is a difference between Random\_Mode\_MAIN\_NI.m (non-intrusive approach: NI-PCE) and Random\_Mode\_MAIN\_Intrusive.m (intrusive approach: I-PCE).
## A. Non intrusive approach 
Three methods are proposed for building the eigenvector PCE. The method is selected at the beginning of Random\_Mode\_MAIN\_NI.m with parameter method:

- method = 1: each element of each eigenvector is identified independently from the other elements; each coefficient identification sample is given in the canonical basis OR in the deterministic modal basis OR in the POD basis.

- method = 2: each eigenvector is identified globally, which means that all the elements of all the samples of one eigenvector are used simultaneously in the eigenvector identification process; each eigenvector identification sample is given in the canonical basis OR in the deterministic modal OR in the POD basis.

The choice of the basis in which the eigenvectors are identified is made at the beginning with parameter choice\_basis

- choice\_basis = 1: canonical basis 

- choice\_basis = 2: modal basis  

- choice\_basis = 3: POD basis, built on the identification samples

The eigenvectors are normalized according to parameter norme\_Psi:

- norme\_Psi = 1: normalization with respect to the mass matrix => modal mass matrix M\_mod=eyes(n\_ddl)

- norme\_Psi = 2: the 1st element of each eigenvector equals 1



But first, the PCE characteristics and the identification and validation samples are determined.


### a. PCE characteristics
In data\_RM\_PC\_NI.m the NI-PCE characteristics are chosen: degree, number of samples to identify the NI-PCE, identification samples, etc.

The identification (and validation) samples are obtained with Random\_Mode\_samples.m, launched from data\_RM.m. 

Therefore, npt\_id sets of uncertain parameters are drawn according to their statistical law.

- Om2\_u\_PC(:,i\_u): PCE squared eigenfrequency vector obtained for the i\_u-th set of uncertain parameters

- Psi\_u\_PC(i\_mod).Psi\_u(:,i\_u): i\_mod-th eigenvector associated with the i\_u-th set of uncertain parameters

    => Psi\_u\_PC(i\_mod).Psi\_u(:,I\_id) : i\_mod-th eigenvector associated with the set of identification data

    => Psi\_u\_PC(i\_mod).Psi\_u(:,I\_val) : i\_mod-th eigenvector associated with the set of validation data


### b. Eigenvector identification
Depending on parameter method, Random\_Mode\_PC\_NI\_BASES.m or Random\_Mode\_PC\_NI\_COORD.m is launched.

Then the squared eigenfrequency, the eigenvectors and the modal masses are identified, with the SVB-ARD method which find a sparse PCE. The sparsity is governed through parameter e\_ARD, which can be tuned in data\_RM.m: the high e\_ARD, the sparser the PCE.

Outputs: 

Squared Eigeneigenfrequencies: PCE\_Om2(i\_mode)

- PCE\_Om(i\_mode).ak: indices of nonzero coefficients
- PCE\_Om(i\_mode).deg\_PC\_sparse: indices of nonzero coefficients


Modal masses: PCE\_M\_mod(i\_mode)

- PCE\_M\_mod (i\_mode).ak: indices of nonzero coefficients
- PCE\_M\_mod (i\_mode).deg\_PC\_sparse: indices of nonzero coefficients


Eigenvectors: PCE\_Psi(i\_mode)

- PCE\_Psi(i\_mode).ak: NI-PCE coefficients of i\_mode-th eigenvector
- PCE\_Psi(i\_mode).Ind\_sparse(i\_Psi\_element).Ind: i\_Psi\_element-th element nonzero coefficients indices of i\_mode-th mode
- PCE\_Psi(i\_mode).deg\_PC\_sparse(i\_Psi\_element).deg: give the polynomials related to the nonzero coefficients

### c. Validation/comparison
The modes and the FRF are calculated with the validation samples set and NI-PCE metamodels.

The results are compared, errors are calculated and plotted

## B. Intrusive approach
### a. PCE characteristics
In data\_RM\_PC\_INT.m the I-PCE characteristics are chosen: degree and maximum number of parameters that are multiplied in each term of the expansion, which is a way to get a sparsity. 


### b. Expectation matrices
The expectations matrices developed in the Appendix of the paper are calculated in Matrices\_moments\_triple\_norm.m, Matrices\_moments\_norm.m and Mat\_mode\_raid\_masse.m


### c. I-PCE coefficients
The coefficients are calculated by solving nonlinear system (A.4) given in Appendix A, using a Newton-Raphson algorithm, which calls nonlinear function Fonction\_nl.m and jacobian Jacobian\_nl.m.

The PCE coefficients are then calculated at the same time:

- PCE\_Psi(i\_mod).lambda: : i\_mod-th eigenvector PCE coefficients in the deterministic modal basis
- PCE\_Om2(i\_mod).ak: i\_mod-th squared eigenfrequency PCE coefficients


### d. Validation/comparison
The modes, the modal masses and the FRF are calculated with the validation samples set and I-PCE metamodels.

The results are compared, errors are calculated and plotted

